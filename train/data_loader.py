"""
Data loading pipeline for RhoFold diffusion training.
Handles .seq, .a3m, and .pdb files for sequence, MSA, and structure data.
"""

import os
import torch
import numpy as np
from torch.utils.data import Dataset, DataLoader
from typing import Dict, List, Tuple, Optional
import pickle
import json
from pathlib import Path

import Bio.PDB
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
import mdtraj as md
from rhofold.utils.alphabet import RNAAlphabet
from rhofold.utils.converter import RNAConverter


class RhoFoldDataset(Dataset):
    """Dataset for loading RhoFold training data."""
    
    def __init__(
        self,
        data_dir: str,
        index_file: str,
        max_seq_len: int = 512,
        msa_depth: int = 128,
        tokenizer=None,
        augment: bool = True,
    ):
        """
        Args:
            data_dir: Directory containing .seq, .a3m, .pdb files
            index_file: Text file with list of sample names (one per line)
            max_seq_len: Maximum sequence length
            msa_depth: Maximum MSA depth
            tokenizer: Tokenizer for sequence encoding
            augment: Whether to apply data augmentation
        """
        self.data_dir = Path(data_dir)
        self.max_seq_len = max_seq_len
        self.msa_depth = msa_depth
        self.augment = augment
        
        # Load sample names
        with open(index_file, 'r') as f:
            self.samples = [line.strip() for line in f if line.strip()]
        
        # Initialize tokenizer
        if tokenizer is None:
            self.alphabet = RNAAlphabet.from_architecture('RNA')
            self.converter = RNAConverter()  # RNAConverter doesn't take alphabet parameter
        else:
            self.alphabet = tokenizer.alphabet
            self.converter = tokenizer.converter
        
        print(f"Loaded {len(self.samples)} samples from {index_file}")
    
    def __len__(self):
        return len(self.samples)
    
    def __getitem__(self, idx):
        sample_name = self.samples[idx]
        
        try:
            # Load sequence
            seq_path = self.data_dir / f"{sample_name}.seq"
            with open(seq_path, 'r') as f:
                sequence = f.read().strip()
            
            # Load MSA
            msa_path = self.data_dir / f"{sample_name}.a3m"
            msa_sequences = self._load_msa(msa_path)
            
            # Load structure
            pdb_path = self.data_dir / f"{sample_name}.pdb"
            structure_data = self._load_structure(pdb_path)
            
            # Process data
            processed = self._process_sample(
                sequence, msa_sequences, structure_data, sample_name
            )
            
            return processed
            
        except Exception as e:
            print(f"Error loading sample {sample_name}: {e}")
            # Return a dummy sample to avoid breaking the batch
            return self._get_dummy_sample()
    
    def _load_msa(self, msa_path: Path) -> List[str]:
        """Load MSA from A3M format."""
        sequences = []
        with open(msa_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                sequences.append(str(record.seq))
        return sequences
    
    def _load_structure(self, pdb_path: Path) -> Dict:
        """Load structure from PDB file."""
        try:
            # Load with mdtraj for coordinate extraction
            traj = md.load(str(pdb_path))
            coords = traj.xyz[0]  # (n_atoms, 3)
            
            # Load with BioPython for detailed structure info
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('structure', str(pdb_path))
            
            # Extract atom coordinates and types
            atom_coords = []
            atom_types = []
            atom_mask = []
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atom_coords.append(atom.get_coord())
                            atom_types.append(self._get_atom_type(atom.get_name()))
                            atom_mask.append(1.0)
            
            atom_coords = np.array(atom_coords)
            atom_types = np.array(atom_types)
            atom_mask = np.array(atom_mask)
            
            return {
                'atom_coords': atom_coords,
                'atom_types': atom_types,
                'atom_mask': atom_mask,
                'traj': traj,
            }
            
        except Exception as e:
            print(f"Error loading structure {pdb_path}: {e}")
            return self._get_dummy_structure()
    
    def _get_atom_type(self, atom_name: str) -> int:
        """Convert atom name to type index."""
        atom_type_map = {
            'C': 0, 'N': 1, 'O': 2, 'P': 3, 'H': 4
        }
        element = atom_name[0] if atom_name else 'C'
        return atom_type_map.get(element, 0)
    
    def _process_sample(self, sequence: str, msa_sequences: List[str], 
                       structure_data: Dict, sample_name: str) -> Dict:
        """Process a single sample into model inputs."""
        
        # Process MSA - tokens should be 2D [msa_depth, seq_len]
        msa_tokens = []
        for seq in msa_sequences[:self.msa_depth]:
            seq_tokens = [self.alphabet.get_idx(s) for s in seq]
            msa_tokens.append(seq_tokens)
        
        # Pad MSA to fixed depth
        while len(msa_tokens) < self.msa_depth:
            msa_tokens.append([0] * len(sequence))  # Pad with zeros
        
        tokens = torch.tensor(msa_tokens, dtype=torch.long)  # [msa_depth, seq_len]
        msa_tokens = tokens  # Same as tokens for now
        
        # Pad sequences to max length
        seq_len = min(tokens.shape[1], self.max_seq_len)
        tokens = tokens[:, :seq_len]  # Keep 2D shape
        msa_tokens = msa_tokens[:, :seq_len]
        
        # Create masks
        msa_mask = torch.ones(msa_tokens.shape[:2], dtype=torch.bool)
        pair_mask = torch.ones((seq_len, seq_len), dtype=torch.bool)
        
        # Process structure data
        atom_coords = torch.tensor(structure_data['atom_coords'], dtype=torch.float32)
        atom_types = torch.tensor(structure_data['atom_types'], dtype=torch.long)
        atom_mask = torch.tensor(structure_data['atom_mask'], dtype=torch.bool)
        
        # Pad structure data to max length
        if len(atom_coords) < self.max_seq_len * 23:  # Assuming 23 atoms per residue
            pad_len = self.max_seq_len * 23 - len(atom_coords)
            atom_coords = torch.cat([
                atom_coords, 
                torch.zeros(pad_len, 3, dtype=torch.float32)
            ])
            atom_types = torch.cat([
                atom_types, 
                torch.zeros(pad_len, dtype=torch.long)
            ])
            atom_mask = torch.cat([
                atom_mask, 
                torch.zeros(pad_len, dtype=torch.bool)
            ])
        
        # Reshape to (batch, seq_len, n_atoms, 3)
        # Ensure we have the right number of atoms per residue
        atoms_per_residue = 23
        total_atoms_needed = seq_len * atoms_per_residue
        
        if len(atom_coords) < total_atoms_needed:
            # Pad with zeros
            pad_len = total_atoms_needed - len(atom_coords)
            atom_coords = torch.cat([
                atom_coords, 
                torch.zeros(pad_len, 3, dtype=torch.float32)
            ])
            atom_types = torch.cat([
                atom_types, 
                torch.zeros(pad_len, dtype=torch.long)
            ])
            atom_mask = torch.cat([
                atom_mask, 
                torch.zeros(pad_len, dtype=torch.bool)
            ])
        elif len(atom_coords) > total_atoms_needed:
            # Truncate
            atom_coords = atom_coords[:total_atoms_needed]
            atom_types = atom_types[:total_atoms_needed]
            atom_mask = atom_mask[:total_atoms_needed]
        
        atom_coords = atom_coords.view(1, seq_len, atoms_per_residue, 3)
        atom_types = atom_types.view(1, seq_len, atoms_per_residue)
        atom_mask = atom_mask.view(1, seq_len, atoms_per_residue)
        
        # Create target data for loss computation
        target = self._create_targets(sequence, structure_data, seq_len)
        
        return {
            'tokens': tokens,
            'msa_tokens': msa_tokens,
            'msa_mask': msa_mask,
            'pair_mask': pair_mask,
            'sequence': sequence,
            'sample_name': sample_name,
            'atom_coords': atom_coords,
            'atom_types': atom_types,
            'atom_mask': atom_mask,
            'target': target,
        }
    
    def _create_targets(self, sequence: str, structure_data: Dict, seq_len: int) -> Dict:
        """Create target data for loss computation."""
        
        # For now, create dummy targets - in real training these would be computed
        # from the structure data
        target = {
            'masked_tokens': torch.full((seq_len,), -100, dtype=torch.long),
            'distance_bins': torch.full((seq_len, seq_len), -100, dtype=torch.long),
            'ss_true': torch.zeros((seq_len, seq_len), dtype=torch.long),
            'atom_coords': structure_data['atom_coords'],
            'atom_mask': structure_data['atom_mask'],
            'pair_mask': torch.ones((seq_len, seq_len), dtype=torch.bool),
            'basepair_list': torch.zeros((0, 2), dtype=torch.long),
            'plddt_true': torch.ones(seq_len, dtype=torch.float32) * 0.8,
            'res_mask': torch.ones(seq_len, dtype=torch.bool),
        }
        
        # Create true frames for FAPE loss
        if 'atom_coords' in structure_data:
            coords = structure_data['atom_coords']
            # Select backbone atoms for frame construction
            # This is a simplified version - in practice you'd select specific atoms
            if len(coords) >= 3:
                selected_coords = coords[:3]  # Take first 3 atoms
                true_frames = self._build_local_frame(selected_coords.reshape(1, 1, 3, 3))
                target['true_frames'] = true_frames
            else:
                target['true_frames'] = {
                    'R': torch.eye(3).unsqueeze(0).unsqueeze(0),
                    't': torch.zeros(1, 1, 3)
                }
        
        return target
    
    def _build_local_frame(self, coords: torch.Tensor) -> Dict:
        """Build local residue frames from coordinates."""
        # Simplified frame construction
        v1 = coords[:, :, 1] - coords[:, :, 0]  # C1' - C4'
        v2 = coords[:, :, 2] - coords[:, :, 0]  # N1/N9 - C4'
        
        # Normalize
        x_axis = v1 / (torch.norm(v1, dim=-1, keepdim=True) + 1e-6)
        
        # Orthogonalize
        proj = (v2 * x_axis).sum(dim=-1, keepdim=True) * x_axis
        y_axis = v2 - proj
        y_axis = y_axis / (torch.norm(y_axis, dim=-1, keepdim=True) + 1e-6)
        
        # Cross product
        z_axis = torch.cross(x_axis, y_axis, dim=-1)
        z_axis = z_axis / (torch.norm(z_axis, dim=-1, keepdim=True) + 1e-6)
        
        # Stack to rotation matrix
        R = torch.stack([x_axis, y_axis, z_axis], dim=-1)
        t = coords[:, :, 0]  # Origin
        
        return {'R': R, 't': t}
    
    def _get_dummy_sample(self) -> Dict:
        """Return a dummy sample for error cases."""
        seq_len = 50
        return {
            'tokens': torch.zeros(seq_len, dtype=torch.long),
            'msa_tokens': torch.zeros(self.msa_depth, seq_len, dtype=torch.long),
            'msa_mask': torch.ones(self.msa_depth, seq_len, dtype=torch.bool),
            'pair_mask': torch.ones(seq_len, seq_len, dtype=torch.bool),
            'sequence': 'A' * seq_len,
            'sample_name': 'dummy',
            'atom_coords': torch.zeros(1, seq_len, 23, 3, dtype=torch.float32),
            'atom_types': torch.zeros(1, seq_len, 23, dtype=torch.long),
            'atom_mask': torch.zeros(1, seq_len, 23, dtype=torch.bool),
            'target': self._get_dummy_target(seq_len),
        }
    
    def _get_dummy_structure(self) -> Dict:
        """Return dummy structure data."""
        return {
            'atom_coords': np.zeros((100, 3), dtype=np.float32),
            'atom_types': np.zeros(100, dtype=np.int64),
            'atom_mask': np.zeros(100, dtype=bool),
        }
    
    def _get_dummy_target(self, seq_len: int) -> Dict:
        """Return dummy target data."""
        return {
            'masked_tokens': torch.full((seq_len,), -100, dtype=torch.long),
            'distance_bins': torch.full((seq_len, seq_len), -100, dtype=torch.long),
            'ss_true': torch.zeros((seq_len, seq_len), dtype=torch.long),
            'atom_coords': torch.zeros(seq_len * 23, 3, dtype=torch.float32),
            'atom_mask': torch.zeros(seq_len * 23, dtype=torch.bool),
            'pair_mask': torch.ones((seq_len, seq_len), dtype=torch.bool),
            'basepair_list': torch.zeros((0, 2), dtype=torch.long),
            'plddt_true': torch.ones(seq_len, dtype=torch.float32) * 0.8,
            'res_mask': torch.ones(seq_len, dtype=torch.bool),
            'true_frames': {
                'R': torch.eye(3).unsqueeze(0).unsqueeze(0),
                't': torch.zeros(1, 1, 3)
            }
        }


def collate_fn(batch):
    """Collate function for DataLoader."""
    # Pad sequences to same length
    max_len = max(item['tokens'].shape[0] for item in batch)
    
    collated = {}
    for key in batch[0].keys():
        if key == 'target':
            # Handle target dict separately
            collated[key] = {}
            for target_key in batch[0]['target'].keys():
                if target_key in ['true_frames']:
                    # Handle frames dict
                    collated[key][target_key] = {
                        'R': torch.stack([item['target'][target_key]['R'] for item in batch]),
                        't': torch.stack([item['target'][target_key]['t'] for item in batch])
                    }
                else:
                    # Pad other tensors
                    tensors = [item['target'][target_key] for item in batch]
                    if tensors[0].dim() == 1:
                        # 1D tensors - check if it's sequence length or flattened
                        if tensors[0].shape[0] == max_len:
                            # Sequence length tensor
                            padded = torch.zeros(len(batch), max_len, dtype=tensors[0].dtype)
                            for i, tensor in enumerate(tensors):
                                padded[i, :len(tensor)] = tensor
                            collated[key][target_key] = padded
                        else:
                            # Flattened tensor (like atom_coords) - just stack
                            collated[key][target_key] = torch.stack(tensors)
                    elif tensors[0].dim() == 2:
                        # 2D tensors - check if it's square (pairwise) or rectangular
                        if tensors[0].shape[0] == tensors[0].shape[1]:
                            # Square matrix (pairwise)
                            padded = torch.zeros(len(batch), max_len, max_len, dtype=tensors[0].dtype)
                            for i, tensor in enumerate(tensors):
                                padded[i, :tensor.shape[0], :tensor.shape[1]] = tensor
                        else:
                            # Rectangular matrix - just stack
                            padded = torch.stack(tensors)
                        collated[key][target_key] = padded
                    else:
                        # Higher dim tensors - stack as is
                        collated[key][target_key] = torch.stack(tensors)
        elif key in ['tokens', 'msa_tokens', 'msa_mask', 'pair_mask']:
            # Pad these tensors
            if key == 'tokens':
                # tokens should be 2D [msa_depth, seq_len]
                if len(batch[0][key].shape) == 2:
                    msa_depth = batch[0][key].shape[0]
                    padded = torch.zeros(len(batch), msa_depth, max_len, dtype=batch[0][key].dtype)
                    for i, item in enumerate(batch):
                        padded[i, :, :item[key].shape[1]] = item[key]
                    collated[key] = padded
                else:
                    # Fallback for 1D tokens - convert to 3D [batch, msa_depth, seq_len]
                    msa_depth = 3  # Default MSA depth
                    padded = torch.zeros(len(batch), msa_depth, max_len, dtype=batch[0][key].dtype)
                    for i, item in enumerate(batch):
                        # Repeat the sequence for all MSA rows
                        for j in range(msa_depth):
                            padded[i, j, :len(item[key])] = item[key]
                    collated[key] = padded
            elif key == 'msa_tokens':
                padded = torch.zeros(len(batch), batch[0][key].shape[0], max_len, dtype=batch[0][key].dtype)
                for i, item in enumerate(batch):
                    padded[i, :, :item[key].shape[1]] = item[key]
                collated[key] = padded
            elif key == 'msa_mask':
                padded = torch.zeros(len(batch), batch[0][key].shape[0], max_len, dtype=batch[0][key].dtype)
                for i, item in enumerate(batch):
                    padded[i, :, :item[key].shape[1]] = item[key]
                collated[key] = padded
            elif key == 'pair_mask':
                padded = torch.zeros(len(batch), max_len, max_len, dtype=batch[0][key].dtype)
                for i, item in enumerate(batch):
                    padded[i, :item[key].shape[0], :item[key].shape[1]] = item[key]
                collated[key] = padded
        elif key in ['atom_coords', 'atom_types', 'atom_mask']:
            # Stack these tensors
            collated[key] = torch.stack([item[key] for item in batch])
        else:
            # For other keys (strings, etc.), just collect them
            collated[key] = [item[key] for item in batch]
    
    return collated


def create_dataloader(
    data_dir: str,
    index_file: str,
    batch_size: int = 4,
    num_workers: int = 4,
    max_seq_len: int = 512,
    msa_depth: int = 128,
    shuffle: bool = True,
) -> DataLoader:
    """Create a DataLoader for RhoFold training."""
    
    dataset = RhoFoldDataset(
        data_dir=data_dir,
        index_file=index_file,
        max_seq_len=max_seq_len,
        msa_depth=msa_depth,
    )
    
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        collate_fn=collate_fn,
        pin_memory=True,
    )
    
    return dataloader


if __name__ == "__main__":
    # Test the data loader
    import tempfile
    import os
    
    # Create dummy data for testing
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create dummy files
        with open(os.path.join(temp_dir, "test.seq"), "w") as f:
            f.write("AUCG" * 10)
        
        with open(os.path.join(temp_dir, "test.a3m"), "w") as f:
            f.write(">seq1\nAUCG" * 10 + "\n>seq2\nAUCG" * 10 + "\n")
        
        # Create dummy PDB
        with open(os.path.join(temp_dir, "test.pdb"), "w") as f:
            f.write("ATOM      1  P   A   1      20.154  16.967  23.862  1.00 20.00           P\n")
            f.write("ATOM      2  O1P A   1      19.123  17.123  24.123  1.00 20.00           O\n")
            f.write("ATOM      3  O2P A   1      21.123  17.123  24.123  1.00 20.00           O\n")
        
        # Create index file
        with open(os.path.join(temp_dir, "index.txt"), "w") as f:
            f.write("test\n")
        
        # Test dataset
        dataset = RhoFoldDataset(temp_dir, os.path.join(temp_dir, "index.txt"))
        print(f"Dataset length: {len(dataset)}")
        
        # Test loading a sample
        sample = dataset[0]
        print(f"Sample keys: {sample.keys()}")
        print(f"Tokens shape: {sample['tokens'].shape}")
        print(f"MSA tokens shape: {sample['msa_tokens'].shape}")
        print(f"Atom coords shape: {sample['atom_coords'].shape}")
        
        # Test dataloader
        dataloader = create_dataloader(temp_dir, os.path.join(temp_dir, "index.txt"), batch_size=1)
        batch = next(iter(dataloader))
        print(f"Batch keys: {batch.keys()}")
        print("Data loader test passed!")
