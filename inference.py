""" """
import logging
import os
import sys

import numpy as np
import torch

from rhofold.data.balstn import BLASTN
from rhofold.rhofold import RhoFold
from rhofold.config import rhofold_config
from rhofold.utils import get_device, save_ss2ct, timing
from rhofold.relax.relax import AmberRelaxation
from rhofold.utils.alphabet import get_features

@torch.no_grad()
def main(config):
    '''
    RhoFold Inference pipeline
    '''

    os.makedirs(config.output_dir, exist_ok=True)

    logger = logging.getLogger('RhoFold Inference')
    logger.setLevel(level=logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')
    file_handler = logging.FileHandler(f'{config.output_dir}/log.txt', mode='w')
    file_handler.setLevel(level=logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    logger.info(f'Constructing RhoFold')
    model = RhoFold(rhofold_config)

    logger.info(f'    loading {config.ckpt}')
    model.load_state_dict(torch.load(config.ckpt, map_location=torch.device('cpu'))['model'])
    model.eval()

    # Input seq, MSA
    logger.info(f"Input_fas {config.input_fas}")

    if config.single_seq_pred:
        config.input_a3m = config.input_fas
        logger.info(f"Input_a3m is None, the modeling will run using single sequence only (input_fas)")

    elif config.input_a3m is None:
        config.input_a3m = f'{config.output_dir}/seq.a3m'
        databases = [f'{config.database_dpath}/rnacentral.fasta', f'{config.database_dpath}/nt']
        blast = BLASTN(binary_dpath=config.binary_dpath, databases=databases)
        blast.query(config.input_fas, f'{config.output_dir}/seq.a3m', logger)
        logger.info(f"Input_a3m {config.input_a3m}")

    else:
        logger.info(f"Input_a3m {config.input_a3m}")

    with timing('RhoFold Inference', logger=logger):

        config.device = get_device(config.device)
        logger.info(f'    Inference using device {config.device}')
        model = model.to(config.device)

        data_dict = get_features(config.input_fas, config.input_a3m)

        # Forward pass
        outputs = model(tokens=data_dict['tokens'].to(config.device),
                       rna_fm_tokens=data_dict['rna_fm_tokens'].to(config.device),
                       seq=data_dict['seq'],
                       )

        output = outputs[-1]

        os.makedirs(config.output_dir, exist_ok=True)

        # Secondary structure, .ct format
        ss_prob_map = torch.sigmoid(output['ss'][0, 0]).data.cpu().numpy()
        ss_file = f'{config.output_dir}/ss.ct'
        save_ss2ct(ss_prob_map, data_dict['seq'], ss_file, threshold=0.5)

        # Dist prob map & Secondary structure prob map, .npz format
        npz_file = f'{config.output_dir}/results.npz'
        np.savez_compressed(npz_file,
                            dist_n = torch.softmax(output['n'].squeeze(0), dim=0).data.cpu().numpy(),
                            dist_p = torch.softmax(output['p'].squeeze(0), dim=0).data.cpu().numpy(),
                            dist_c = torch.softmax(output['c4_'].squeeze(0), dim=0).data.cpu().numpy(),
                            ss_prob_map = ss_prob_map,
                            plddt = output['plddt'][0].data.cpu().numpy(),
                            )

        # Save the prediction
        unrelaxed_model = f'{config.output_dir}/unrelaxed_model.pdb'

        # The last cords prediction
        node_cords_pred = output['cord_tns_pred'][-1].squeeze(0)
        model.structure_module.converter.export_pdb_file(data_dict['seq'],
                                                         node_cords_pred.data.cpu().numpy(),
                                                         path=unrelaxed_model, chain_id=None,
                                                         confidence=output['plddt'][0].data.cpu().numpy(),
                                                         logger=logger)

    # Amber relaxation
    if config.relax_steps is not None:
        relax_steps = int(config.relax_steps)
        if relax_steps > 0:
            with timing(f'Amber Relaxation : {relax_steps} iterations', logger=logger):
                amber_relax = AmberRelaxation(max_iterations=relax_steps, logger=logger)
                relaxed_model = f'{config.output_dir}/relaxed_{relax_steps}_model.pdb'
                amber_relax.process(unrelaxed_model, relaxed_model)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--device", help="Default cpu. If GPUs are available, you can set --device cuda:<GPU_index> for faster prediction.", default=None)
    parser.add_argument("--ckpt", help="Path to the pretrained model, default ./pretrained/model_20221010_params.pt", default='./pretrained/rhofold_pretrained_params.pt')
    parser.add_argument("--input_fas", help="Path to the input fasta file. Valid nucleic acids in RNA sequence: A, U, G, C", required=True)
    parser.add_argument("--input_a3m", help="Path to the input msa file. Default None."
                                            "If --input_a3m is not given (set to None), MSA will be generated automatically. ", default=None)
    parser.add_argument("--output_dir", help="Path to the output dir. "
                                             "3D prediction is saved in .pdb format. "
                                             "Distogram prediction is saved in .npz format. "
                                             "Secondary structure prediction is save in .ct format. ", required=True)
    parser.add_argument("--relax_steps", help="Num of steps for structure refinement, default 1000.", default = 1000)
    parser.add_argument("--single_seq_pred", help="Default False. If --single_seq_pred is set to True, "
                                                       "the modeling will run using single sequence only (input_fas)", default=False)
    parser.add_argument("--database_dpath", help="Path to the pretrained model, default ./database", default='./database')
    parser.add_argument("--binary_dpath", help="Path to the pretrained model, default ./rhofold/data/bin", default='./rhofold/data/bin')

    args = parser.parse_args()
    main(args)
