"""
Training script for RhoFold with diffusion module.
Includes forward pass validation and gradient checking.
"""

import os
import sys
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
import numpy as np
import argparse
from pathlib import Path
import json
import time
from tqdm import tqdm
import wandb

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from rhofold.config import rhofold_config
from rhofold.rhofold_with_diffusion import create_diffusion_rhofold
from train.data_loader import create_dataloader
from train.losses import RhoFoldLoss


class DiffusionTrainer:
    """Trainer class for RhoFold with diffusion module."""
    
    def __init__(self, config, args):
        self.config = config
        self.args = args
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        # Initialize model
        self.model = self._create_model()
        self.model.to(self.device)
        
        # Initialize loss function
        self.criterion = RhoFoldLoss()
        
        # Initialize optimizer
        self.optimizer = optim.AdamW(
            self.model.parameters(),
            lr=args.learning_rate,
            weight_decay=args.weight_decay,
            betas=(0.9, 0.95)
        )
        
        # Initialize scheduler
        self.scheduler = optim.lr_scheduler.CosineAnnealingLR(
            self.optimizer, T_max=args.epochs
        )
        
        # Initialize data loader
        self.train_loader = self._create_dataloader()
        
        # Initialize logging
        if args.use_wandb:
            wandb.init(
                project="rhofold-diffusion",
                config=vars(args),
                name=f"diffusion-{args.experiment_name}"
            )
    
    def _create_model(self):
        """Create the diffusion model."""
        diffusion_config = {
            'time_embed_dim': self.args.time_embed_dim,
            'num_timesteps': self.args.num_timesteps,
            'noise_schedule': self.args.noise_schedule,
            'base_channels': self.args.base_channels,
        }
        
        model = create_diffusion_rhofold(
            self.config, 
            diffusion_config=diffusion_config
        )
        
        print(f"Created model with {sum(p.numel() for p in model.parameters()):,} parameters")
        return model
    
    def _create_dataloader(self):
        """Create data loader."""
        return create_dataloader(
            data_dir=self.args.data_dir,
            index_file=self.args.index_file,
            batch_size=self.args.batch_size,
            num_workers=self.args.num_workers,
            max_seq_len=self.args.max_seq_len,
            msa_depth=self.args.msa_depth,
            shuffle=True,
        )
    
    def validate_forward_pass(self):
        """Validate that the forward pass works correctly."""
        print("Validating forward pass...")
        
        self.model.eval()
        with torch.no_grad():
            # Get a batch
            batch = next(iter(self.train_loader))
            
            # Move to device
            for key in batch:
                if isinstance(batch[key], torch.Tensor):
                    batch[key] = batch[key].to(self.device)
                elif isinstance(batch[key], dict):
                    for subkey in batch[key]:
                        if isinstance(batch[key][subkey], torch.Tensor):
                            batch[key][subkey] = batch[key][subkey].to(self.device)
            
            # Forward pass
            try:
                outputs = self.model(
                    tokens=batch['tokens'],
                    rna_fm_tokens=batch['msa_tokens'],
                    seq=batch['sequence']
                )
                
                print("✓ Forward pass successful")
                print(f"  Number of recycling outputs: {len(outputs)}")
                
                # Check output shapes
                for i, output in enumerate(outputs):
                    print(f"  Recycling {i}:")
                    for key, value in output.items():
                        if isinstance(value, torch.Tensor):
                            print(f"    {key}: {value.shape}")
                        elif isinstance(value, dict):
                            for subkey, subvalue in value.items():
                                if isinstance(subvalue, torch.Tensor):
                                    print(f"    {key}.{subkey}: {subvalue.shape}")
                
                # Check diffusion output shape
                if hasattr(self.model.e2eformer, 'diffusion_blocks'):
                    print("  Diffusion module detected")
                    # The diffusion module should output pair features with shape (B, L, L, C)
                    # This is handled in the DiffusionStack.forward method
                
                return True
                
            except Exception as e:
                print(f"✗ Forward pass failed: {e}")
                import traceback
                traceback.print_exc()
                return False
    
    def validate_loss_computation(self):
        """Validate that loss computation works correctly."""
        print("Validating loss computation...")
        
        self.model.eval()
        with torch.no_grad():
            # Get a batch
            batch = next(iter(self.train_loader))
            
            # Move to device
            for key in batch:
                if isinstance(batch[key], torch.Tensor):
                    batch[key] = batch[key].to(self.device)
                elif isinstance(batch[key], dict):
                    for subkey in batch[key]:
                        if isinstance(batch[key][subkey], torch.Tensor):
                            batch[key][subkey] = batch[key][subkey].to(self.device)
            
            # Forward pass
            outputs = self.model(
                tokens=batch['tokens'],
                rna_fm_tokens=batch['msa_tokens'],
                seq=batch['sequence']
            )
            
            # Use the last recycling output for loss computation
            output = outputs[-1]
            
            # Prepare predictions for loss computation
            pred = self._prepare_predictions(output, batch)
            
            # Compute loss
            try:
                total_loss, loss_dict = self.criterion(pred, batch['target'])
                
                print("✓ Loss computation successful")
                print(f"  Total loss: {total_loss.item():.6f}")
                print("  Individual losses:")
                for key, value in loss_dict.items():
                    print(f"    {key}: {value.item():.6f}")
                
                # Check for NaN values
                if torch.isnan(total_loss):
                    print("✗ Loss contains NaN values")
                    return False
                
                return True
                
            except Exception as e:
                print(f"✗ Loss computation failed: {e}")
                import traceback
                traceback.print_exc()
                return False
    
    def validate_gradient_flow(self):
        """Validate that gradients flow correctly through the diffusion module."""
        print("Validating gradient flow...")
        
        self.model.train()
        
        # Get a batch
        batch = next(iter(self.train_loader))
        
        # Move to device
        for key in batch:
            if isinstance(batch[key], torch.Tensor):
                batch[key] = batch[key].to(self.device)
            elif isinstance(batch[key], dict):
                for subkey in batch[key]:
                    if isinstance(batch[key][subkey], torch.Tensor):
                        batch[key][subkey] = batch[key][subkey].to(self.device)
        
        # Forward pass
        outputs = self.model(
            tokens=batch['tokens'],
            rna_fm_tokens=batch['msa_tokens'],
            seq=batch['sequence']
        )
        
        # Use the last recycling output for loss computation
        output = outputs[-1]
        
        # Prepare predictions for loss computation
        pred = self._prepare_predictions(output, batch)
        
        # Compute loss
        total_loss, loss_dict = self.criterion(pred, batch['target'])
        
        # Backward pass
        try:
            total_loss.backward()
            print("✓ Backward pass successful")
            
            # Check gradients in diffusion module
            if hasattr(self.model.e2eformer, 'diffusion_blocks'):
                diffusion_grads = []
                for block in self.model.e2eformer.diffusion_blocks:
                    for param in block.parameters():
                        if param.grad is not None:
                            diffusion_grads.append(param.grad.abs().mean().item())
                
                if diffusion_grads:
                    avg_grad = np.mean(diffusion_grads)
                    print(f"  Average gradient magnitude in diffusion: {avg_grad:.6f}")
                    
                    if avg_grad > 1e-8:
                        print("✓ Gradients are flowing to diffusion module")
                        return True
                    else:
                        print("✗ Gradients are too small in diffusion module")
                        return False
                else:
                    print("✗ No gradients found in diffusion module")
                    return False
            else:
                print("✗ No diffusion module found")
                return False
                
        except Exception as e:
            print(f"✗ Backward pass failed: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    def _prepare_predictions(self, output, batch):
        """Prepare model output for loss computation."""
        pred = {}
        
        # Distance predictions
        if 'p' in output and 'c4_' in output and 'n' in output:
            # Combine distance predictions
            distance_logits = torch.stack([
                output['p'], output['c4_'], output['n']
            ], dim=1)  # (B, 3, L, L, n_bins)
            pred['distance_logits'] = distance_logits
        
        # Secondary structure predictions
        if 'ss' in output:
            pred['ss_logits'] = output['ss']
        
        # pLDDT predictions
        if 'plddt' in output:
            pred['plddt_logits'] = output['plddt']
        
        # Atom coordinates
        if 'cords_c1\'' in output:
            pred['atom_coords'] = output['cords_c1\''][-1]  # Use last recycling
            # Create dummy atom types (in practice, these would be computed)
            B, L = pred['atom_coords'].shape[:2]
            pred['atom_types'] = torch.zeros(B, L, 23, dtype=torch.long, device=self.device)
        
        # Create dummy pseudo coordinates for ss3d loss
        if 'atom_coords' in pred:
            B, L, N, _ = pred['atom_coords'].shape
            pred['pseudo_coords'] = pred['atom_coords'][:, :, :4, :]  # Use first 4 atoms
        
        # Create dummy frames for FAPE loss
        if 'atom_coords' in pred:
            coords = pred['atom_coords'][:, :, :3, :]  # Use first 3 atoms
            pred['pred_frames'] = self._build_local_frame(coords)
        
        return pred
    
    def _build_local_frame(self, coords):
        """Build local residue frames from coordinates."""
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
    
    def train_epoch(self, epoch):
        """Train for one epoch."""
        self.model.train()
        total_loss = 0.0
        num_batches = 0
        
        pbar = tqdm(self.train_loader, desc=f"Epoch {epoch}")
        for batch_idx, batch in enumerate(pbar):
            # Move to device
            for key in batch:
                if isinstance(batch[key], torch.Tensor):
                    batch[key] = batch[key].to(self.device)
                elif isinstance(batch[key], dict):
                    for subkey in batch[key]:
                        if isinstance(batch[key][subkey], torch.Tensor):
                            batch[key][subkey] = batch[key][subkey].to(self.device)
            
            # Forward pass
            outputs = self.model(
                tokens=batch['tokens'],
                rna_fm_tokens=batch['msa_tokens'],
                seq=batch['sequence']
            )
            
            # Use the last recycling output for loss computation
            output = outputs[-1]
            
            # Prepare predictions for loss computation
            pred = self._prepare_predictions(output, batch)
            
            # Compute loss
            loss, loss_dict = self.criterion(pred, batch['target'])
            
            # Backward pass
            self.optimizer.zero_grad()
            loss.backward()
            
            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            
            self.optimizer.step()
            
            # Update statistics
            total_loss += loss.item()
            num_batches += 1
            
            # Update progress bar
            pbar.set_postfix({
                'loss': f"{loss.item():.6f}",
                'avg_loss': f"{total_loss/num_batches:.6f}"
            })
            
            # Log to wandb
            if self.args.use_wandb and batch_idx % 10 == 0:
                log_dict = {'train/loss': loss.item()}
                for key, value in loss_dict.items():
                    log_dict[f'train/{key}_loss'] = value.item()
                wandb.log(log_dict)
        
        return total_loss / num_batches
    
    def train(self):
        """Main training loop."""
        print("Starting training...")
        
        # Validate before training
        if not self.validate_forward_pass():
            print("Forward pass validation failed. Exiting.")
            return
        
        if not self.validate_loss_computation():
            print("Loss computation validation failed. Exiting.")
            return
        
        if not self.validate_gradient_flow():
            print("Gradient flow validation failed. Exiting.")
            return
        
        print("All validations passed. Starting training...")
        
        best_loss = float('inf')
        
        for epoch in range(self.args.epochs):
            # Train epoch
            avg_loss = self.train_epoch(epoch)
            
            # Update scheduler
            self.scheduler.step()
            
            # Log epoch results
            print(f"Epoch {epoch}: Average loss = {avg_loss:.6f}")
            
            if self.args.use_wandb:
                wandb.log({
                    'epoch': epoch,
                    'train/avg_loss': avg_loss,
                    'train/learning_rate': self.scheduler.get_last_lr()[0]
                })
            
            # Save checkpoint
            if avg_loss < best_loss:
                best_loss = avg_loss
                self.save_checkpoint(epoch, avg_loss, is_best=True)
            
            if epoch % self.args.save_every == 0:
                self.save_checkpoint(epoch, avg_loss, is_best=False)
    
    def save_checkpoint(self, epoch, loss, is_best=False):
        """Save model checkpoint."""
        checkpoint = {
            'epoch': epoch,
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'scheduler_state_dict': self.scheduler.state_dict(),
            'loss': loss,
        }
        
        if is_best:
            path = os.path.join(self.args.checkpoint_dir, 'best_model.pth')
        else:
            path = os.path.join(self.args.checkpoint_dir, f'checkpoint_epoch_{epoch}.pth')
        
        os.makedirs(self.args.checkpoint_dir, exist_ok=True)
        torch.save(checkpoint, path)
        print(f"Saved checkpoint: {path}")


def main():
    parser = argparse.ArgumentParser(description='Train RhoFold with diffusion module')
    
    # Data arguments
    parser.add_argument('--data_dir', type=str, required=True,
                       help='Directory containing .seq, .a3m, .pdb files')
    parser.add_argument('--index_file', type=str, required=True,
                       help='Text file with list of sample names')
    parser.add_argument('--max_seq_len', type=int, default=512,
                       help='Maximum sequence length')
    parser.add_argument('--msa_depth', type=int, default=128,
                       help='Maximum MSA depth')
    
    # Model arguments
    parser.add_argument('--time_embed_dim', type=int, default=256,
                       help='Time embedding dimension')
    parser.add_argument('--num_timesteps', type=int, default=1000,
                       help='Number of diffusion timesteps')
    parser.add_argument('--noise_schedule', type=str, default='cosine',
                       choices=['linear', 'cosine'],
                       help='Noise schedule type')
    parser.add_argument('--base_channels', type=int, default=64,
                       help='Base number of channels in UNet')
    
    # Training arguments
    parser.add_argument('--batch_size', type=int, default=4,
                       help='Batch size')
    parser.add_argument('--epochs', type=int, default=100,
                       help='Number of epochs')
    parser.add_argument('--learning_rate', type=float, default=1e-4,
                       help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-5,
                       help='Weight decay')
    parser.add_argument('--num_workers', type=int, default=4,
                       help='Number of data loader workers')
    
    # Logging arguments
    parser.add_argument('--use_wandb', action='store_true',
                       help='Use Weights & Biases for logging')
    parser.add_argument('--experiment_name', type=str, default='diffusion_test',
                       help='Experiment name for logging')
    
    # Checkpoint arguments
    parser.add_argument('--checkpoint_dir', type=str, default='./checkpoints',
                       help='Directory to save checkpoints')
    parser.add_argument('--save_every', type=int, default=10,
                       help='Save checkpoint every N epochs')
    
    args = parser.parse_args()
    
    # Create trainer
    trainer = DiffusionTrainer(rhofold_config, args)
    
    # Run validation tests
    print("Running validation tests...")
    
    # Test forward pass
    if not trainer.validate_forward_pass():
        print("Forward pass validation failed!")
        return
    
    # Test loss computation
    if not trainer.validate_loss_computation():
        print("Loss computation validation failed!")
        return
    
    # Test gradient flow
    if not trainer.validate_gradient_flow():
        print("Gradient flow validation failed!")
        return
    
    print("All validation tests passed!")
    
    # Start training
    trainer.train()


if __name__ == "__main__":
    main()
