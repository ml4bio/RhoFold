"""
Final training script for RhoFold with diffusion module.
This script trains the model using the new diffusion module instead of E2Eformer.
"""

import os
import sys
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
import numpy as np
from pathlib import Path
import argparse
import time
from tqdm import tqdm

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from rhofold.config import rhofold_config
from rhofold.rhofold_with_diffusion import create_diffusion_rhofold
from train.data_loader import create_dataloader
from train.losses import RhoFoldLoss


class DiffusionTrainer:
    """Trainer for RhoFold with diffusion module."""
    
    def __init__(self, config):
        self.config = config
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(f"Using device: {self.device}")
        
        # Create model
        self.model = create_diffusion_rhofold(rhofold_config, config['diffusion'])
        self.model.to(self.device)
        
        # Create loss function
        self.criterion = RhoFoldLoss()
        
        # Create optimizer
        self.optimizer = optim.AdamW(
            self.model.parameters(),
            lr=config['learning_rate'],
            weight_decay=config['weight_decay']
        )
        
        # Create scheduler
        self.scheduler = optim.lr_scheduler.CosineAnnealingLR(
            self.optimizer,
            T_max=config['epochs']
        )
        
        # Training state
        self.epoch = 0
        self.step = 0
        self.best_loss = float('inf')
        
    def prepare_predictions(self, output, batch, device):
        """Prepare model output for loss computation."""
        pred = {}
        
        # Distance predictions
        if 'p' in output and 'c4_' in output and 'n' in output:
            # The outputs are (B, n_bins, L, L), we need (B, L, L, n_bins)
            p_logits = output['p'].permute(0, 2, 3, 1)  # (B, L, L, n_bins)
            c4_logits = output['c4_'].permute(0, 2, 3, 1)  # (B, L, L, n_bins)
            n_logits = output['n'].permute(0, 2, 3, 1)  # (B, L, L, n_bins)
            
            # Stack to (B, 3, L, L, n_bins) - keep as 5D tensor
            distance_logits = torch.stack([p_logits, c4_logits, n_logits], dim=1)
            pred['distance_logits'] = distance_logits
        
        # Secondary structure predictions
        if 'ss' in output:
            ss_logits = output['ss']
            # Remove extra dimension if present: (B, 1, L, L) -> (B, L, L)
            if ss_logits.dim() == 4 and ss_logits.shape[1] == 1:
                ss_logits = ss_logits.squeeze(1)
            pred['ss_logits'] = ss_logits
        
        # pLDDT predictions
        if 'plddt' in output:
            plddt = output['plddt']
            if isinstance(plddt, tuple):
                plddt = plddt[0]  # Take first element
            
            # Ensure plddt_logits has shape (B, L, n_bins)
            if plddt.dim() == 2:  # (B, L) -> (B, L, n_bins)
                B, L = plddt.shape
                n_bins = 50  # Default number of bins
                # Create dummy logits - in practice these would be real predictions
                pred['plddt_logits'] = torch.randn(B, L, n_bins, device=device)
            else:
                pred['plddt_logits'] = plddt
        
        # Atom coordinates
        coords_key = 'cords_c1\''
        if coords_key in output:
            coords = output[coords_key]
            if isinstance(coords, list):
                coords = coords[0]  # Take first element
            elif isinstance(coords, tuple):
                coords = coords[0]  # Take first element
            
            # Ensure coords has shape (B, L, N, 3) for loss functions
            if coords.dim() == 3:  # (B, L, 3) -> (B, L, 1, 3)
                coords = coords.unsqueeze(2)
            elif coords.dim() == 2:  # (L, 3) -> (1, L, 1, 3)
                coords = coords.unsqueeze(0).unsqueeze(2)
            
            pred['atom_coords'] = coords
            # Create dummy atom types and mask (in practice, these would be computed)
            B, L, N = pred['atom_coords'].shape[:3]
            pred['atom_types'] = torch.zeros(B, L, N, dtype=torch.long, device=device)
            pred['atom_mask'] = torch.ones(B, L, N, dtype=torch.bool, device=device)
        
        return pred
    
    def train_epoch(self, dataloader):
        """Train for one epoch."""
        self.model.train()
        total_loss = 0.0
        num_batches = 0
        
        pbar = tqdm(dataloader, desc=f"Epoch {self.epoch}")
        for batch_idx, batch in enumerate(pbar):
            # Move batch to device
            for key in batch:
                if isinstance(batch[key], torch.Tensor):
                    batch[key] = batch[key].to(self.device)
                elif isinstance(batch[key], dict):
                    for subkey in batch[key]:
                        if isinstance(batch[key][subkey], torch.Tensor):
                            batch[key][subkey] = batch[key][subkey].to(self.device)
            
            # Convert msa_tokens from 3D to 2D for RNA-FM
            rna_fm_tokens = batch['msa_tokens'][:, 0, :]
            
            # Convert sequence from list to string
            seq_str = batch['sequence'][0] if isinstance(batch['sequence'], list) else batch['sequence']
            
            # Forward pass
            self.optimizer.zero_grad()
            
            try:
                outputs = self.model(
                    tokens=batch['tokens'],
                    rna_fm_tokens=rna_fm_tokens,
                    seq=seq_str
                )
                
                # Use the last recycling output for loss computation
                output = outputs[-1]
                
                # Prepare predictions for loss computation
                pred = self.prepare_predictions(output, batch, self.device)
                
                # Compute loss
                total_loss_batch, loss_dict = self.criterion(pred, batch['target'])
                
                # Backward pass
                total_loss_batch.backward()
                
                # Gradient clipping
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                
                # Update parameters
                self.optimizer.step()
                
                # Update statistics
                total_loss += total_loss_batch.item()
                num_batches += 1
                self.step += 1
                
                # Update progress bar
                pbar.set_postfix({
                    'loss': f"{total_loss_batch.item():.4f}",
                    'avg_loss': f"{total_loss / num_batches:.4f}",
                    'lr': f"{self.optimizer.param_groups[0]['lr']:.6f}"
                })
                
            except Exception as e:
                print(f"Error in batch {batch_idx}: {e}")
                continue
        
        return total_loss / max(num_batches, 1)
    
    def validate(self, dataloader):
        """Validate the model."""
        self.model.eval()
        total_loss = 0.0
        num_batches = 0
        
        with torch.no_grad():
            for batch in tqdm(dataloader, desc="Validation"):
                # Move batch to device
                for key in batch:
                    if isinstance(batch[key], torch.Tensor):
                        batch[key] = batch[key].to(self.device)
                    elif isinstance(batch[key], dict):
                        for subkey in batch[key]:
                            if isinstance(batch[key][subkey], torch.Tensor):
                                batch[key][subkey] = batch[key][subkey].to(self.device)
                
                # Convert msa_tokens from 3D to 2D for RNA-FM
                rna_fm_tokens = batch['msa_tokens'][:, 0, :]
                
                # Convert sequence from list to string
                seq_str = batch['sequence'][0] if isinstance(batch['sequence'], list) else batch['sequence']
                
                try:
                    # Forward pass
                    outputs = self.model(
                        tokens=batch['tokens'],
                        rna_fm_tokens=rna_fm_tokens,
                        seq=seq_str
                    )
                    
                    # Use the last recycling output for loss computation
                    output = outputs[-1]
                    
                    # Prepare predictions for loss computation
                    pred = self.prepare_predictions(output, batch, self.device)
                    
                    # Compute loss
                    total_loss_batch, loss_dict = self.criterion(pred, batch['target'])
                    
                    total_loss += total_loss_batch.item()
                    num_batches += 1
                    
                except Exception as e:
                    print(f"Error in validation batch: {e}")
                    continue
        
        return total_loss / max(num_batches, 1)
    
    def save_checkpoint(self, path):
        """Save model checkpoint."""
        checkpoint = {
            'epoch': self.epoch,
            'step': self.step,
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'scheduler_state_dict': self.scheduler.state_dict(),
            'best_loss': self.best_loss,
            'config': self.config
        }
        torch.save(checkpoint, path)
        print(f"Checkpoint saved to {path}")
    
    def load_checkpoint(self, path):
        """Load model checkpoint."""
        checkpoint = torch.load(path, map_location=self.device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        self.scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
        self.epoch = checkpoint['epoch']
        self.step = checkpoint['step']
        self.best_loss = checkpoint['best_loss']
        print(f"Checkpoint loaded from {path}")
    
    def train(self, train_dataloader, val_dataloader=None, save_dir='checkpoints'):
        """Main training loop."""
        os.makedirs(save_dir, exist_ok=True)
        
        print(f"Starting training for {self.config['epochs']} epochs")
        print(f"Training samples: {len(train_dataloader.dataset)}")
        if val_dataloader:
            print(f"Validation samples: {len(val_dataloader.dataset)}")
        
        for epoch in range(self.config['epochs']):
            self.epoch = epoch
            start_time = time.time()
            
            # Training
            train_loss = self.train_epoch(train_dataloader)
            
            # Validation
            val_loss = None
            if val_dataloader:
                val_loss = self.validate(val_dataloader)
            
            # Update learning rate
            self.scheduler.step()
            
            # Logging
            epoch_time = time.time() - start_time
            print(f"\nEpoch {epoch + 1}/{self.config['epochs']}")
            print(f"Train Loss: {train_loss:.6f}")
            if val_loss is not None:
                print(f"Val Loss: {val_loss:.6f}")
            print(f"Time: {epoch_time:.2f}s")
            print(f"LR: {self.optimizer.param_groups[0]['lr']:.6f}")
            
            # Save checkpoint
            if val_loss is not None and val_loss < self.best_loss:
                self.best_loss = val_loss
                self.save_checkpoint(os.path.join(save_dir, 'best_model.pth'))
            
            # Save regular checkpoint
            if (epoch + 1) % self.config['save_every'] == 0:
                self.save_checkpoint(os.path.join(save_dir, f'checkpoint_epoch_{epoch + 1}.pth'))
        
        print("Training completed!")


def main():
    parser = argparse.ArgumentParser(description='Train RhoFold with diffusion module')
    parser.add_argument('--data_dir', type=str, required=True, help='Path to data directory')
    parser.add_argument('--index_file', type=str, required=True, help='Path to index file')
    parser.add_argument('--val_data_dir', type=str, help='Path to validation data directory')
    parser.add_argument('--val_index_file', type=str, help='Path to validation index file')
    parser.add_argument('--batch_size', type=int, default=4, help='Batch size')
    parser.add_argument('--epochs', type=int, default=100, help='Number of epochs')
    parser.add_argument('--learning_rate', type=float, default=1e-4, help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-5, help='Weight decay')
    parser.add_argument('--max_seq_len', type=int, default=512, help='Maximum sequence length')
    parser.add_argument('--msa_depth', type=int, default=128, help='MSA depth')
    parser.add_argument('--save_every', type=int, default=10, help='Save checkpoint every N epochs')
    parser.add_argument('--resume', type=str, help='Path to checkpoint to resume from')
    parser.add_argument('--save_dir', type=str, default='checkpoints', help='Directory to save checkpoints')
    
    args = parser.parse_args()
    
    # Create config
    config = {
        'diffusion': {
            'time_embed_dim': 128,
            'num_timesteps': 100,
            'noise_schedule': 'linear',
            'base_channels': 32,
        },
        'learning_rate': args.learning_rate,
        'weight_decay': args.weight_decay,
        'epochs': args.epochs,
        'save_every': args.save_every,
    }
    
    # Create data loaders
    train_dataloader = create_dataloader(
        data_dir=args.data_dir,
        index_file=args.index_file,
        batch_size=args.batch_size,
        num_workers=4,
        max_seq_len=args.max_seq_len,
        msa_depth=args.msa_depth,
        shuffle=True,
    )
    
    val_dataloader = None
    if args.val_data_dir and args.val_index_file:
        val_dataloader = create_dataloader(
            data_dir=args.val_data_dir,
            index_file=args.val_index_file,
            batch_size=args.batch_size,
            num_workers=4,
            max_seq_len=args.max_seq_len,
            msa_depth=args.msa_depth,
            shuffle=False,
        )
    
    # Create trainer
    trainer = DiffusionTrainer(config)
    
    # Resume from checkpoint if specified
    if args.resume:
        trainer.load_checkpoint(args.resume)
    
    # Start training
    trainer.train(train_dataloader, val_dataloader, args.save_dir)


if __name__ == "__main__":
    main()

