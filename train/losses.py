# rhofold_losses.py
import torch
import torch.nn as nn
import torch.nn.functional as F

# ---------- utilities ----------
def mask_mean(x, mask, dim=None, eps=1e-8):
    """Masked mean. mask: (..., L, L) or (..., L) -> expand as needed."""
    if mask is None:
        return x.mean(dim=dim)
    # ensure mask dtype float
    m = mask.to(x.dtype)
    if dim is None:
        return (x * m).sum() / (m.sum() + eps)
    else:
        return (x * m).sum(dim=dim) / (m.sum(dim=dim) + eps)


# ---------- 1D: MLM loss ----------
def L_mlm(logits, target_tokens, mask):
    """
    logits: (B, L, vocab) - predicted logits for MLM
    target_tokens: (B, L) - integer token ids (with -100 for masked / ignored)
    mask: (B, L) boolean where positions to compute loss are True (the 5% masked)
    returns: scalar loss (CE averaged over masked positions)
    """
    # use cross_entropy with ignore_index
    logits_flat = logits.view(-1, logits.size(-1))
    targets_flat = target_tokens.view(-1)
    loss = F.cross_entropy(logits_flat, targets_flat, reduction='mean', ignore_index=-100)
    return loss


# ---------- 2D: distance loss ----------
def L_dis(distance_logits, true_dists_bins, pair_mask):
    """
    distance_logits: (B, 3, L, L, n_bins) or (B, L, L, n_bins) logits over distance bins for each pair
    true_dists_bins: (B, L, L) integer bin indices (0..n_bins-1), use -100 for ignore
    pair_mask: (B, L, L) boolean mask for valid pair entries
    """
    if distance_logits.dim() == 5:
        # Handle (B, 3, L, L, n_bins) - process each atom type separately
        B, n_types, L, _, n_bins = distance_logits.shape
        total_loss = 0.0
        for i in range(n_types):
            logits = distance_logits[:, i].view(-1, n_bins)  # (B*L*L, n_bins)
            targets = true_dists_bins.view(-1)  # (B*L*L,)
            loss = F.cross_entropy(logits, targets, ignore_index=-100, reduction='none')  # (B*L*L,)
            loss = loss.view(B, L, L)
            if pair_mask is not None:
                loss = loss * pair_mask.to(loss.dtype)
                denom = pair_mask.sum(dim=(1,2)).clamp(min=1.0)
                total_loss += loss.sum(dim=(1,2)).mean()  # average across batch
            else:
                total_loss += loss.mean()
        return total_loss / n_types  # average across atom types
    else:
        # Handle (B, L, L, n_bins) - original case
        B, L, _, n_bins = distance_logits.shape
        logits = distance_logits.view(-1, n_bins)  # (B*L*L, n_bins)
        targets = true_dists_bins.view(-1)  # (B*L*L,)
        loss = F.cross_entropy(logits, targets, ignore_index=-100, reduction='none')  # (B*L*L,)
        loss = loss.view(B, L, L)
        if pair_mask is not None:
            loss = loss * pair_mask.to(loss.dtype)
            denom = pair_mask.sum(dim=(1,2)).clamp(min=1.0)
            return loss.sum(dim=(1,2)).mean()  # average across batch
        else:
            return loss.mean()


# ---------- 2D: secondary structure (binary) ----------
def L_ss(ss_logits, ss_true, pair_mask):
    """
    ss_logits: (B, L, L) logits (or (B,L,L,2) if two-class logits)
    ss_true: (B, L, L) 0/1
    pair_mask: (B, L, L)
    """
    # support both single-logit or two-class
    if ss_logits.dim() == 4 and ss_logits.size(-1) == 2:
        logits = ss_logits.view(-1, 2)
        targets = ss_true.view(-1)
        loss = F.cross_entropy(logits, targets, ignore_index=-100, reduction='none').view_as(ss_true)
    else:
        # use BCEWithLogits
        loss = F.binary_cross_entropy_with_logits(ss_logits, ss_true.to(ss_logits.dtype), reduction='none')
    if pair_mask is not None:
        loss = loss * pair_mask.to(loss.dtype)
        return loss.sum() / (pair_mask.sum() + 1e-8)
    else:
        return loss.mean()


# ---------- 3D: FAPE loss (Frame Aligned Point Error) ----------
def build_local_frame(coords, atom_mask=None, eps=1e-6):
    """
    Build local residue frames from 3 atoms: (C4', C1', N1/N9).
    coords: (B, L, 3_atoms, 3)  -- selected backbone atoms per residue
    atom_mask: (B, L, 3_atoms) optional mask for missing atoms
    returns dict { 'R': (B,L,3,3), 't': (B,L,3) }
    """
    v1 = coords[:,:,1] - coords[:,:,0]   # (B,L,3) C1' - C4'
    v2 = coords[:,:,2] - coords[:,:,0]   # (B,L,3) N1/N9 - C4'

    # normalize
    x_axis = v1 / (torch.norm(v1, dim=-1, keepdim=True) + eps)

    # orthogonalize v2 against x_axis
    proj = (v2 * x_axis).sum(dim=-1, keepdim=True) * x_axis
    y_axis = v2 - proj
    y_axis = y_axis / (torch.norm(y_axis, dim=-1, keepdim=True) + eps)

    # z axis = x × y
    z_axis = torch.cross(x_axis, y_axis, dim=-1)
    z_axis = z_axis / (torch.norm(z_axis, dim=-1, keepdim=True) + eps)

    # stack
    R = torch.stack([x_axis, y_axis, z_axis], dim=-1)  # (B,L,3,3)
    t = coords[:,:,0]  # use C4' as origin (B,L,3)

    return {'R': R, 't': t}

def L_fape(pred_frames, pred_atom_coords, true_frames, true_atom_coords, atom_mask, clamp_distance=10.0, eps=1e-8):
    """
    Frame Aligned Point Error (FAPE)
    pred_frames: dict {'R': (B,L,3,3), 't': (B,L,3)}
    true_frames: dict {'R': (B,L,3,3), 't': (B,L,3)}
    pred_atom_coords, true_atom_coords: (B,L,N,3)
    atom_mask: (B,L,N) boolean
    clamp_distance: clip distance (like AlphaFold, prevents exploding loss)
    returns: scalar loss
    """
    R_p, t_p = pred_frames['R'], pred_frames['t']  # (B,L,3,3), (B,L,3)
    R_t, t_t = true_frames['R'], true_frames['t']

    # local coords: transform atoms into each frame
    x_pred = pred_atom_coords - t_p.unsqueeze(2)   # (B,L,N,3)
    x_true = true_atom_coords - t_t.unsqueeze(2)

    # apply rotations: R^T * (x - t)
    x_pred_local = torch.einsum('blij,blnj->blni', R_p.transpose(-2,-1), x_pred)
    x_true_local = torch.einsum('blij,blnj->blni', R_t.transpose(-2,-1), x_true)

    # distance error
    d = torch.norm(x_pred_local - x_true_local, dim=-1)  # (B,L,N)

    if clamp_distance is not None:
        d = torch.clamp(d, max=clamp_distance)

    if atom_mask is not None:
        m = atom_mask.to(d.dtype)
        loss = (d * m).sum() / (m.sum() + eps)
    else:
        loss = d.mean()

    return loss


# ---------- 3D: clash loss ----------
def L_clash(atom_coords, atom_types, atom_mask, vdw_radii=None, clash_tol=0.4):
    """
    Penalize atom pairs that are closer than sum of vdw radii minus tolerance.
    atom_coords: (B, L, N, 3)
    atom_types: (B, L, N) integer mapping to element types (e.g., map to 'C','N','O','P' indices)
    atom_mask: (B, L, N)
    vdw_radii: dict or tensor mapping type_idx -> radius (float)
    clash_tol: margin (if dist < rsum - clash_tol => penalize)
    Implementation details: O(N^2) per-residue naive; we compute pairwise across all atoms per structure.
    """
    B, L, N, _ = atom_coords.shape
    device = atom_coords.device
    if vdw_radii is None:
        # default small radii for typical atoms (mapping 0:C,1:N,2:O,3:P,4:H)
        default = torch.tensor([1.7, 1.55, 1.52, 1.80, 1.2], device=device)
    else:
        default = torch.tensor([vdw_radii[i] for i in range(len(vdw_radii))], device=device)

    # flatten per model
    coords = atom_coords.view(B, L*N, 3)
    mask = atom_mask.view(B, L*N)
    types = atom_types.view(B, L*N)

    # pairwise distances (B, M, M)
    d = torch.cdist(coords, coords)  # (B, M, M)
    # rsum matrix
    radii = default[types]  # (B, M)
    rsum = radii.unsqueeze(2) + radii.unsqueeze(1)  # (B, M, M)
    # compute violation: where d < rsum - clash_tol
    thresh = (rsum - clash_tol)
    # exclude self pairs by large mask
    large = 1e6 * torch.eye(d.size(-1), device=device).unsqueeze(0)
    d = d + large

    violation = F.relu(thresh - d)  # positive if violating
    # mask out invalid atoms
    m = mask.unsqueeze(2) * mask.unsqueeze(1)  # (B, M, M)
    violation = violation * m.to(violation.dtype)
    # sum and average
    loss = violation.sum() / (m.sum() + 1e-8)
    return loss


# ---------- 3D: Lss3d (pseudo-atom base pairing constraint) ----------
def L_ss3d(pseudo_coords, basepair_list, dm_stats=None, tau=0.5):
    """
    pseudo_coords: (B, L, 4, 3) four pseudo atoms per residue in local frame (or global coords)
    basepair_list: list or tensor of shape (B, Npairs, 2) indices (m,n) that are base-paired in this sample
    dm_stats: if provided, statistical dm for pair types, dict keyed by base pair type -> (4,4) distances
              else dm computed from ground truth coords (should be passed in target)
    tau: tolerance
    Implementation follows equation in paper: max(pred - tau - dm, 0) summed over pairs & pseudo atoms.
    """
    B, L, _, _ = pseudo_coords.shape
    device = pseudo_coords.device
    loss_total = 0.0
    total_pairs = 0
    for b in range(B):
        pairs = basepair_list[b]  # (Np, 2) maybe padded with -1
        valid = pairs[:,0] >= 0
        if valid.sum() == 0:
            continue
        pairs = pairs[valid]
        total_pairs += pairs.shape[0]
        for (m, n) in pairs:
            # predicted distances between four pseudo atoms: (4,4)
            pm = pseudo_coords[b, m]  # (4,3)
            pn = pseudo_coords[b, n]
            d_pred = torch.cdist(pm.unsqueeze(0), pn.unsqueeze(0))[0]  # (4,4)
            # get dm (standard) either from stats keyed by pair type or passed array
            if dm_stats is None:
                # fallback: use d_pred as target (zero loss) -- but in practice pass true dm
                dm = d_pred.detach()
            else:
                # assuming dm_stats[b][m,n] or keyed by base pair type
                dm = dm_stats[b][m,n]  # (4,4), precomputed
            # compute hinge: max(d_pred - tau - dm, 0)
            hinge = F.relu(d_pred - tau - dm)
            loss_total += hinge.sum()
    if total_pairs == 0:
        return torch.tensor(0.0, device=device)
    return loss_total / float(total_pairs)


# ---------- pLDDT loss ----------
def L_pLDDT(plddt_logits, plddt_true, mask=None, n_bins=50):
    """
    plddt_logits: (B, L, n_bins) logits predicting binned LDDT (bins width 0.02)
    plddt_true: (B, L) float in [0,1] or actual LDDT in 0..100 depending on preprocessing; map accordingly
    mask: (B, L)
    """
    # map true value to bin index
    # assume plddt_true in [0,1]; if in [0,100] convert: plddt_true/100
    if plddt_true.max() > 2.0:
        p_true = (plddt_true / 100.0).clamp(0.0, 1.0)
    else:
        p_true = plddt_true.clamp(0.0, 1.0)
    bin_idx = torch.clamp((p_true * (n_bins - 1)).long(), 0, n_bins - 1)  # (B, L)
    logits = plddt_logits.view(-1, n_bins)
    targets = bin_idx.view(-1)
    loss = F.cross_entropy(logits, targets, reduction='none', ignore_index=-100)
    loss = loss.view_as(bin_idx).to(logits.dtype)
    if mask is not None:
        loss = (loss * mask.to(loss.dtype)).sum() / (mask.sum() + 1e-8)
    else:
        loss = loss.mean()
    return loss


# ---------- full combined loss wrapper ----------
class RhoFoldLoss(nn.Module):
    def __init__(self, weights=None):
        super().__init__()
        # default weights as in paper
        if weights is None:
            self.weights = {
                'mlm': 1.0,
                'dis': 0.3,
                'ss': 0.1,
                'clash': 0.03,
                'fape': 2.0,
                'ss3d': 0.1,
                'plddt': 0.01
            }
        else:
            self.weights = weights

    def forward(self, pred, target):
        """
        pred: dict with keys:
            - 'mlm_logits' (B,L,V)
            - 'distance_logits' (B,L,L,n_bins)
            - 'ss_logits' (B,L,L) or (B,L,L,2)
            - 'atom_coords' (B,L,N,3)
            - 'pred_frames' {'R':(B,L,3,3),'t':(B,L,3)}
            - 'plddt_logits' (B,L,n_bins)
            - 'pseudo_coords' (B,L,4,3)
            - 'atom_types' (B,L,N)
        target: dict with keys:
            - 'masked_tokens' (B,L) with -100 for ignore
            - 'distance_bins' (B,L,L) with -100 ignore
            - 'ss_true' (B,L,L)
            - 'atom_coords' (B,L,N,3)
            - 'true_frames' dict like pred_frames
            - 'atom_mask' (B,L,N)
            - 'pair_mask' (B,L,L)
            - 'basepair_list' (B, Npairs, 2) or (B,L,L) adjacency
            - 'plddt_true' (B,L)
            - 'dm_stats' optional for ss3d
        """
        losses = {}
        # MLM
        if 'mlm_logits' in pred and 'masked_tokens' in target:
            losses['mlm'] = L_mlm(pred['mlm_logits'], target['masked_tokens'], None)
        else:
            # Use device from pred tensors since RhoFoldLoss has no parameters
            device = next(iter(pred.values())).device if pred else torch.device('cpu')
            losses['mlm'] = torch.tensor(0.0, device=device)

        losses['dis'] = L_dis(pred['distance_logits'], target['distance_bins'], target.get('pair_mask', None))
        losses['ss'] = L_ss(pred['ss_logits'], target['ss_true'], target.get('pair_mask', None))

        losses['clash'] = L_clash(pred['atom_coords'], pred['atom_types'], pred.get('atom_mask', None))

        # FAPE: require frames; fallback to coordinate RMS if frames missing
        if 'pred_frames' in pred and 'true_frames' in target:
            # Ensure target atom_coords has the same shape as pred atom_coords
            true_coords = target['atom_coords']
            if true_coords.dim() == 5:  # (B, 1, L, N, 3) -> (B, L, N, 3)
                true_coords = true_coords.squeeze(1)
            elif true_coords.dim() == 4 and true_coords.shape[1] == 1:  # (B, 1, L, 3) -> (B, L, 1, 3)
                true_coords = true_coords.unsqueeze(2)
            
            # Ensure shapes match
            if true_coords.shape != pred['atom_coords'].shape:
                # Pad or truncate to match
                B, L, N, _ = pred['atom_coords'].shape
                if true_coords.dim() == 3:  # (B, L, 3) -> (B, L, 1, 3)
                    true_coords = true_coords.unsqueeze(2)
                elif true_coords.dim() == 4:
                    if true_coords.shape[2] < N:
                        # Pad with zeros
                        pad_size = N - true_coords.shape[2]
                        true_coords = torch.cat([true_coords, torch.zeros(B, L, pad_size, 3, device=true_coords.device)], dim=2)
                    elif true_coords.shape[2] > N:
                        # Truncate
                        true_coords = true_coords[:, :, :N, :]
            
            # Ensure true_frames has correct shape
            true_frames = target['true_frames']
            if true_frames['R'].dim() == 4 and true_frames['R'].shape[0] == 1 and true_frames['R'].shape[1] == 1:
                # (1, 1, 3, 3) -> (B, L, 3, 3)
                B, L = pred['atom_coords'].shape[:2]
                true_frames = {
                    'R': true_frames['R'].expand(B, L, 3, 3),
                    't': true_frames['t'].expand(B, L, 3)
                }
            
            losses['fape'] = L_fape(pred['pred_frames'], pred['atom_coords'], true_frames, true_coords, pred.get('atom_mask', None))
        else:
            # fallback: RMSD per atom masked
            true_coords = target['atom_coords']
            if true_coords.dim() == 5:  # (B, 1, L, N, 3) -> (B, L, N, 3)
                true_coords = true_coords.squeeze(1)
            elif true_coords.dim() == 4 and true_coords.shape[1] == 1:  # (B, 1, L, 3) -> (B, L, 1, 3)
                true_coords = true_coords.unsqueeze(2)
            
            # Ensure shapes match
            if true_coords.shape != pred['atom_coords'].shape:
                B, L, N, _ = pred['atom_coords'].shape
                if true_coords.dim() == 3:  # (B, L, 3) -> (B, L, 1, 3)
                    true_coords = true_coords.unsqueeze(2)
                elif true_coords.dim() == 4:
                    if true_coords.shape[2] < N:
                        pad_size = N - true_coords.shape[2]
                        true_coords = torch.cat([true_coords, torch.zeros(B, L, pad_size, 3, device=true_coords.device)], dim=2)
                    elif true_coords.shape[2] > N:
                        true_coords = true_coords[:, :, :N, :]
            
            diff = (pred['atom_coords'] - true_coords).pow(2).sum(dim=-1)
            m = pred.get('atom_mask', None)
            if m is not None:
                losses['fape'] = (diff * m.to(diff.dtype)).sum() / (m.sum() + 1e-8)
            else:
                losses['fape'] = diff.mean()

        losses['ss3d'] = L_ss3d(pred.get('pseudo_coords', None), target.get('basepair_list', None), target.get('dm_stats', None))

        losses['plddt'] = L_pLDDT(pred['plddt_logits'], target['plddt_true'], mask=target.get('res_mask', None))

        # combine with weights
        total = 0.0
        for k, w in self.weights.items():
            total = total + w * losses[k]
        return total, losses