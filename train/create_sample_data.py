"""
Create sample data for testing RhoFold diffusion training.
This script generates dummy .seq, .a3m, and .pdb files for testing.
"""

import os
import argparse
import numpy as np
from pathlib import Path


def create_rna_sequence(length, base_probs=None):
    """Create a random RNA sequence."""
    if base_probs is None:
        base_probs = [0.25, 0.25, 0.25, 0.25]  # Equal probability for A, U, C, G
    
    bases = ['A', 'U', 'C', 'G']
    return ''.join(np.random.choice(bases, size=length, p=base_probs))


def create_msa_sequences(sequence, num_sequences=5, mutation_rate=0.1):
    """Create MSA sequences by adding mutations to the base sequence."""
    sequences = [sequence]  # First sequence is the original
    
    for _ in range(num_sequences - 1):
        mutated = list(sequence)
        for i in range(len(mutated)):
            if np.random.random() < mutation_rate:
                # Randomly mutate to another base
                bases = ['A', 'U', 'C', 'G']
                bases.remove(mutated[i])  # Remove current base
                mutated[i] = np.random.choice(bases)
        sequences.append(''.join(mutated))
    
    return sequences


def create_pdb_structure(sequence, output_path):
    """Create a simple PDB structure for the RNA sequence."""
    with open(output_path, 'w') as f:
        f.write("HEADER    GENERATED RNA STRUCTURE\n")
        
        atom_id = 1
        residue_id = 1
        
        for i, base in enumerate(sequence):
            # Create a simple helix-like structure
            x = i * 3.4  # 3.4 Angstrom spacing
            y = 0.0
            z = 0.0
            
            # Add some variation
            y += np.sin(i * 0.5) * 2.0
            z += np.cos(i * 0.3) * 1.5
            
            # Create backbone atoms
            atoms = [
                ('P', [x, y, z]),
                ('O1P', [x-1, y+1, z]),
                ('O2P', [x+1, y+1, z]),
                ('O5\'', [x, y-1, z]),
                ('C5\'', [x, y-2, z]),
                ('C4\'', [x, y-3, z]),
                ('O4\'', [x, y-4, z]),
                ('C1\'', [x, y-5, z]),
                ('N1', [x, y-6, z]),
                ('C2', [x, y-7, z]),
                ('N3', [x, y-8, z]),
                ('C4', [x, y-9, z]),
                ('C5', [x, y-10, z]),
                ('C6', [x, y-11, z]),
                ('N6', [x, y-12, z]),
                ('N7', [x, y-13, z]),
                ('C8', [x, y-14, z]),
                ('N9', [x, y-15, z]),
                ('C2\'', [x, y-16, z]),
                ('C3\'', [x, y-17, z]),
                ('O3\'', [x, y-18, z]),
                ('H1', [x, y-19, z]),
                ('H2', [x, y-20, z]),
            ]
            
            for atom_name, (ax, ay, az) in atoms:
                f.write(f"ATOM  {atom_id:5d}  {atom_name:4s} A{residue_id:4d}    {ax:8.3f}{ay:8.3f}{az:8.3f}  1.00 20.00           {atom_name[0]}\n")
                atom_id += 1
            
            residue_id += 1
        
        f.write("END\n")


def create_sample_data(output_dir, num_samples=10, seq_length=100, msa_depth=5):
    """Create sample data for testing."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create index file
    index_file = output_dir / "index.txt"
    
    print(f"Creating {num_samples} sample sequences...")
    
    with open(index_file, 'w') as f:
        for i in range(num_samples):
            sample_name = f"sample_{i:03d}"
            f.write(f"{sample_name}\n")
            
            # Create sequence
            sequence = create_rna_sequence(seq_length)
            
            # Create .seq file
            seq_file = output_dir / f"{sample_name}.seq"
            with open(seq_file, 'w') as sf:
                sf.write(sequence)
            
            # Create .a3m file
            msa_file = output_dir / f"{sample_name}.a3m"
            msa_sequences = create_msa_sequences(sequence, msa_depth)
            with open(msa_file, 'w') as mf:
                for j, seq in enumerate(msa_sequences):
                    mf.write(f">{sample_name}_seq_{j}\n{seq}\n")
            
            # Create .pdb file
            pdb_file = output_dir / f"{sample_name}.pdb"
            create_pdb_structure(sequence, pdb_file)
            
            if (i + 1) % 10 == 0:
                print(f"Created {i + 1}/{num_samples} samples")
    
    print(f"✓ Created {num_samples} samples in {output_dir}")
    print(f"✓ Index file: {index_file}")
    
    # Create a summary
    summary_file = output_dir / "summary.txt"
    with open(summary_file, 'w') as f:
        f.write("Sample Data Summary\n")
        f.write("==================\n")
        f.write(f"Number of samples: {num_samples}\n")
        f.write(f"Sequence length: {seq_length}\n")
        f.write(f"MSA depth: {msa_depth}\n")
        f.write(f"Files per sample: .seq, .a3m, .pdb\n")
        f.write(f"Index file: index.txt\n")
    
    print(f"✓ Summary file: {summary_file}")


def main():
    parser = argparse.ArgumentParser(description='Create sample data for RhoFold training')
    parser.add_argument('--output_dir', type=str, default='./sample_data',
                       help='Output directory for sample data')
    parser.add_argument('--num_samples', type=int, default=10,
                       help='Number of samples to create')
    parser.add_argument('--seq_length', type=int, default=100,
                       help='Length of each sequence')
    parser.add_argument('--msa_depth', type=int, default=5,
                       help='Number of sequences in MSA')
    
    args = parser.parse_args()
    
    print("RhoFold Sample Data Generator")
    print("============================")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of samples: {args.num_samples}")
    print(f"Sequence length: {args.seq_length}")
    print(f"MSA depth: {args.msa_depth}")
    print()
    
    create_sample_data(
        output_dir=args.output_dir,
        num_samples=args.num_samples,
        seq_length=args.seq_length,
        msa_depth=args.msa_depth
    )
    
    print("\nSample data creation completed!")
    print(f"You can now test the training pipeline with:")
    print(f"  python train/validate_setup.py")
    print(f"  python train/test_forward.py")


if __name__ == "__main__":
    main()
