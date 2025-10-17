# Modified RhoFold with diffusion module integration
# This shows how to replace E2EformerStack with DiffusionStack

import torch
import torch.nn as nn

from rhofold.model.embedders import MSAEmbedder, RecyclingEmbedder
# from rhofold.model.e2eformer import E2EformerStack  # Original transformer
from rhofold.model.diffusion import DiffusionStack  # New diffusion module
from rhofold.model.structure_module import StructureModule
from rhofold.model.heads import DistHead, SSHead, pLDDTHead
from rhofold.utils.tensor_utils import add
from rhofold.utils import exists


class RhoFoldWithDiffusion(nn.Module):
    """RhoFold network with diffusion module replacing E2Eformer."""

    def __init__(self, config, use_diffusion=True, diffusion_config=None):
        """Constructor function."""
        super().__init__()

        self.config = config
        self.use_diffusion = use_diffusion

        self.msa_embedder = MSAEmbedder(
            **config.model.msa_embedder,
        )
        
        # Choose between transformer and diffusion
        if use_diffusion:
            # Use diffusion module
            if diffusion_config is None:
                diffusion_config = {
                    'time_embed_dim': 256,
                    'num_timesteps': 1000,
                    'noise_schedule': 'cosine',
                    'base_channels': 64,
                }
            
            self.e2eformer = DiffusionStack(
                **config.model.e2eformer_stack,
                diffusion_config=diffusion_config,
            )
        else:
            # Use original transformer (would need to import E2EformerStack)
            from rhofold.model.e2eformer import E2EformerStack
            self.e2eformer = E2EformerStack(
                **config.model.e2eformer_stack,
            )
        
        self.structure_module = StructureModule(
            **config.model.structure_module,
        )
        self.recycle_embnet = RecyclingEmbedder(
            **config.model.recycling_embedder,
        )
        self.dist_head = DistHead(
            **config.model.heads.dist,
        )
        self.ss_head = SSHead(
            **config.model.heads.ss,
        )
        self.plddt_head = pLDDTHead(
            **config.model.heads.plddt,
        )

    def forward_cords(self, tokens, single_fea, pair_fea, seq):
        output = self.structure_module.forward(seq, tokens, { "single": single_fea, "pair": pair_fea } )
        output['plddt'] = self.plddt_head(output['single'][-1])
        return output

    def forward_heads(self, pair_fea):
        output = {}
        output['ss'] = self.ss_head(pair_fea.float())
        output['p'], output['c4_'], output['n'] = self.dist_head(pair_fea.float())
        return output

    def forward_one_cycle(self, tokens, rna_fm_tokens, recycling_inputs, seq, valid_msa=None):
        '''
        Args:
            tokens: [bs, seq_len, c_z]
            rna_fm_tokens: [bs, seq_len, c_z]
            valid_msa: Boolean tensor indicating which samples have valid MSA [bs]
        '''

        device = tokens.device

        msa_tokens_pert = tokens[:, :self.config.globals.msa_depth]

        msa_fea, pair_fea = self.msa_embedder.forward(tokens=msa_tokens_pert,
                                                      rna_fm_tokens=rna_fm_tokens,
                                                      is_BKL=True)
        
        # Handle invalid MSA samples
        if valid_msa is not None:
            # For samples without valid MSA, zero out MSA embeddings
            # and rely solely on RNA-FM sequence embeddings
            invalid_msa_mask = ~valid_msa
            if invalid_msa_mask.any():
                # Zero out MSA features for invalid samples
                msa_fea[invalid_msa_mask] = 0.0
                # Keep pair features as they come from RNA-FM

        if exists(self.recycle_embnet) and exists(recycling_inputs):
            msa_fea_up, pair_fea_up = self.recycle_embnet(recycling_inputs['single_fea'],
                                                          recycling_inputs['pair_fea'],
                                                          recycling_inputs["cords_c1'"])
            msa_fea[..., 0, :, :] += msa_fea_up
            pair_fea = add(pair_fea, pair_fea_up, inplace=False)

        # Use diffusion or transformer for processing
        msa_fea, pair_fea, single_fea = self.e2eformer(
            m=msa_fea,
            z=pair_fea,
            msa_mask=torch.ones(msa_fea.shape[:3]).to(device),
            pair_mask=torch.ones(pair_fea.shape[:3]).to(device),
            chunk_size=None,
        )

        output = self.forward_cords(tokens, single_fea, pair_fea, seq)
        output.update(self.forward_heads(pair_fea))

        recycling_outputs = {
            'single_fea': msa_fea[..., 0, :, :].detach(),
            'pair_fea': pair_fea.detach(),
            "cords_c1'": output["cords_c1'"][-1].detach(),
        }

        return output, recycling_outputs

    def forward(self,
                tokens,
                rna_fm_tokens,
                seq,
                valid_msa=None,
                **kwargs):

        """Perform the forward pass.

        Args:
            tokens: MSA tokens [B, MSA_depth, L]
            rna_fm_tokens: RNA-FM tokens [B, L]
            seq: Sequence string
            valid_msa: Boolean tensor indicating which samples have valid MSA [B]

        Returns:
        """

        recycling_inputs = None

        outputs = []
        for _r in range(self.config.model.recycling_embedder.recycles):
            output, recycling_inputs = \
                self.forward_one_cycle(tokens, rna_fm_tokens, recycling_inputs, seq, valid_msa)
            outputs.append(output)

        return outputs


def create_diffusion_rhofold(config, diffusion_config=None):
    """Create RhoFold with diffusion module."""
    return RhoFoldWithDiffusion(config, use_diffusion=True, diffusion_config=diffusion_config)


def create_original_rhofold(config):
    """Create original RhoFold with transformer."""
    return RhoFoldWithDiffusion(config, use_diffusion=False)


# Example usage
if __name__ == "__main__":
    from rhofold.config import rhofold_config
    
    # Create diffusion configuration
    diffusion_config = {
        'time_embed_dim': 256,
        'num_timesteps': 1000,
        'noise_schedule': 'cosine',
        'base_channels': 64,
    }
    
    # Create RhoFold with diffusion
    model = create_diffusion_rhofold(rhofold_config, diffusion_config)
    
    print("RhoFold with diffusion module created successfully!")
    print(f"Model uses diffusion: {model.use_diffusion}")
    print(f"Number of parameters: {sum(p.numel() for p in model.parameters()):,}")
