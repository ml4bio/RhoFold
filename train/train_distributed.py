"""
Distributed training script for RhoFold with diffusion module.
Supports multi-GPU training on 4x2080Ti setup.
"""

import os
import sys
import torch
import torch.nn as nn
import torch.distributed as dist
import torch.multiprocessing as mp
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data.distributed import DistributedSampler
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


def setup_distributed(rank, world_size, backend='nccl'):
    """Initialize distributed training."""
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = '12355'
    
    # Initialize the process group
    dist.init_process_group(backend, rank=rank, world_size=world_size)
    
    # Set device
    torch.cuda.set_device(rank)


def cleanup_distributed():
    """Clean up distributed training."""
    dist.destroy_process_group()


class DistributedDiffusionTrainer:
    """Distributed trainer class for RhoFold with diffusion module."""
    
    def __init__(self, config, args, rank, world_size):
        self.config = config
        self.args = args
        self.rank = rank
        self.world_size = world_size
        self.device = torch.device(f'cuda:{rank}')
        
        # Initialize model
        self.model = self._create_model()
        self.model.to(self.device)
        
        # Wrap with DDP
        self.model = DDP(
            self.model,
            device_ids=[rank],
            output_device=rank,
            find_unused_parameters=True
        )
        
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
        
        # Initialize logging (only on rank 0)
        if args.use_wandb and rank == 0:
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
        
        if self.rank == 0:
            print(f"Created model with {sum(p.numel() for p in model.parameters()):,} parameters")
        return model
    
    def _create_dataloader(self):
        """Create distributed data loader."""
        # Create dataset
        from train.data_loader import RhoFoldDataset
        dataset = RhoFoldDataset(
            data_dir=self.args.data_dir,
            index_file=self.args.index_file,
            max_seq_len=self.args.max_seq_len,
            msa_depth=self.args.msa_depth,
        )
        
        # Create distributed sampler
        sampler = DistributedSampler(
            dataset,
            num_replicas=self.world_size,
            rank=self.rank,
            shuffle=True
        )
        
        # Create data loader
        from train.data_loader import collate_fn
        dataloader = torch.utils.data.DataLoader(
            dataset,
            batch_size=self.args.batch_size,
            sampler=sampler,
            num_workers=self.args.num_workers,
            collate_fn=collate_fn,
            pin_memory=True,
            drop_last=True,  # Important for distributed training
        )
        
        return dataloader
    
    def train_epoch(self, epoch):
        """Train for one epoch."""
        self.model.train()
        self.train_loader.sampler.set_epoch(epoch)  # Important for distributed training
        
        total_loss = 0.0
        num_batches = 0
        
        if self.rank == 0:
            pbar = tqdm(self.train_loader, desc=f"Epoch {epoch}")
        else:
            pbar = self.train_loader
        
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
            
            # Update progress bar (only on rank 0)
            if self.rank == 0:
                pbar.set_postfix({
                    'loss': f"{loss.item():.6f}",
                    'avg_loss': f"{total_loss/num_batches:.6f}"
                })
            
            # Log to wandb (only on rank 0)
            if self.args.use_wandb and self.rank == 0 and batch_idx % 10 == 0:
                log_dict = {'train/loss': loss.item()}
                for key, value in loss_dict.items():
                    log_dict[f'train/{key}_loss'] = value.item()
                wandb.log(log_dict)
        
        return total_loss / num_batches
    
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
    
    def train(self):
        """Main training loop."""
        if self.rank == 0:
            print("Starting distributed training...")
        
        best_loss = float('inf')
        
        for epoch in range(self.args.epochs):
            # Train epoch
            avg_loss = self.train_epoch(epoch)
            
            # Update scheduler
            self.scheduler.step()
            
            # Log epoch results (only on rank 0)
            if self.rank == 0:
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
            'model_state_dict': self.model.module.state_dict(),  # Remove DDP wrapper
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


def train_worker(rank, world_size, config, args):
    """Worker function for distributed training."""
    # Setup distributed
    setup_distributed(rank, world_size)
    
    # Create trainer
    trainer = DistributedDiffusionTrainer(config, args, rank, world_size)
    
    # Train
    trainer.train()
    
    # Cleanup
    cleanup_distributed()


def main():
    parser = argparse.ArgumentParser(description='Distributed training for RhoFold with diffusion module')
    
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
                       help='Batch size per GPU')
    parser.add_argument('--epochs', type=int, default=100,
                       help='Number of epochs')
    parser.add_argument('--learning_rate', type=float, default=1e-4,
                       help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-5,
                       help='Weight decay')
    parser.add_argument('--num_workers', type=int, default=4,
                       help='Number of data loader workers per GPU')
    
    # Distributed training arguments
    parser.add_argument('--world_size', type=int, default=4,
                       help='Number of GPUs to use')
    parser.add_argument('--backend', type=str, default='nccl',
                       help='Distributed backend')
    
    # Logging arguments
    parser.add_argument('--use_wandb', action='store_true',
                       help='Use Weights & Biases for logging')
    parser.add_argument('--experiment_name', type=str, default='diffusion_distributed',
                       help='Experiment name for logging')
    
    # Checkpoint arguments
    parser.add_argument('--checkpoint_dir', type=str, default='./checkpoints',
                       help='Directory to save checkpoints')
    parser.add_argument('--save_every', type=int, default=10,
                       help='Save checkpoint every N epochs')
    
    args = parser.parse_args()
    
    # Check if CUDA is available
    if not torch.cuda.is_available():
        print("CUDA is not available. Please run on a GPU-enabled machine.")
        return
    
    # Check number of available GPUs
    available_gpus = torch.cuda.device_count()
    if available_gpus < args.world_size:
        print(f"Warning: Only {available_gpus} GPUs available, but {args.world_size} requested.")
        args.world_size = available_gpus
    
    print(f"Starting distributed training on {args.world_size} GPUs")
    print(f"Batch size per GPU: {args.batch_size}")
    print(f"Total effective batch size: {args.batch_size * args.world_size}")
    
    # Start distributed training
    mp.spawn(
        train_worker,
        args=(args.world_size, rhofold_config, args),
        nprocs=args.world_size,
        join=True
    )


if __name__ == "__main__":
    main()
