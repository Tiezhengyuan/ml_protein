'''
unconditional design
'''
import re
import os
import time
import pickle
import torch
import numpy as np
import random
import glob
from omegaconf import OmegaConf
import py3Dmol

from rfdiffusion.util import writepdb_multi, writepdb
from rfdiffusion.chemical import seq2chars

class RunRfd:
            
    def __init__(self, sampler):
        self.sampler = sampler
        self.device = torch.cuda.get_device_name(torch.cuda.current_device()) \
            if torch.cuda.is_available() else "CPU"

    def run(self, overwrite=False):
        '''
        run RFDiffusion
        '''
        for i, out_prefix, traj_prefix in self._iter(overwrite):
            self.make_deterministic(i)
            start_time = time.time()

            print(f"Making design {out_prefix}")
            self.train()

            # pX0 last step
            # Now don't output sidechains
            writepdb(
                f"{out_prefix}.pdb",
                self.denoised_xyz_stack[0, :, :4],
                self.final_seq,
                self.sampler.binderlen,
                chain_idx=self.sampler.chain_idx,
                bfacts=self.bfacts,
            )
            # run metadata
            self.to_trb(out_prefix, start_time)

            if traj_prefix:
                # trajectory pdbs
                writepdb_multi(
                    f"{traj_prefix}_Xt-1_traj.pdb",
                    self.denoised_xyz_stack,
                    self.bfacts,
                    self.final_seq.squeeze(),
                    use_hydrogens=False,
                    backbone_only=False,
                    chain_ids=self.sampler.chain_idx,
                )
                writepdb_multi(
                    f"{traj_prefix}_pX0_traj.pdb",
                    self.px0_xyz_stack,
                    self.bfacts,
                    self.final_seq.squeeze(),
                    use_hydrogens=False,
                    backbone_only=False,
                    chain_ids=self.sampler.chain_idx,
                )
            print(f"Finished design in {(time.time() - start_time)/60:.2f} minutes")

    def _iter(self, overwrite:bool):
        '''
        set prefix for exports
        '''
        start = self.sampler._conf.inference.design_startnum
        end = start + self.sampler.inf_conf.num_designs
        for i_des in range(start, end):
            # prefix for <output>
            out_prefix = f"{self.sampler.inf_conf.output_prefix}_{i_des}"
            if overwrite is False and self.sampler.inf_conf.cautious \
                and os.path.exists(out_prefix + ".pdb"):
                print(f"(cautious mode) Skipping this design because {out_prefix}.pdb already exists.")
            else:
                # Save outputs
                os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

                # prefix for /traj
                traj_prefix = None
                if self.sampler.inf_conf.write_trajectory:
                    # trajectory pdbs
                    traj_prefix = (
                        os.path.dirname(out_prefix) + "/traj/" + os.path.basename(out_prefix)
                    )
                    os.makedirs(os.path.dirname(traj_prefix), exist_ok=True)
                yield i_des, out_prefix, traj_prefix

    def make_deterministic(self, seed=0):
        if self.sampler._conf.inference.deterministic:
            torch.manual_seed(seed)
            np.random.seed(seed)
            random.seed(seed)

    def train(self):
        '''
        train RFDiffusion model
        '''
        denoised_xyz_stack = []
        # px0[t+1] is provided as a template input to the model at time t
        # px0: (L,14,3) The model's prediction of x0.
        px0_xyz_stack = []
        seq_stack = []
        # plddt: (L, 1) Predicted lDDT of x0.
        plddt_stack = []
        
        # seq_init (L,22) The initialized sequence used in updating the sequence.
        x_init, seq_init = self.sampler.sample_init()
        # x_t (L,14,3) The residue positions at the beginning of this timestep
        # x_t_1: (L,14,3) The updated positions of the next step.
        x_t = torch.clone(x_init)
        # seq_t (L,22) The sequence at the beginning of this timestep
        # seq_t_1: (L) The sequence to the next step (== seq_init)
        seq_t = torch.clone(seq_init)
        # Loop over number of reverse diffusion time steps.
        start = int(self.sampler.t_step_input)
        end = self.sampler.inf_conf.final_step - 1
        # t The timestep that has just been predicted
        for t in range(start, end, -1):
            # Generate the next pos that the model should be supplied at timestep t-1
            px0, x_t, seq_t, plddt = self.sampler.sample_step(
                t=t,
                x_t=x_t,
                seq_init=seq_t,
                final_step=self.sampler.inf_conf.final_step
            )
            px0_xyz_stack.append(px0)
            denoised_xyz_stack.append(x_t)
            seq_stack.append(seq_t)
            # remove singleton leading dimension
            plddt_stack.append(plddt[0])

        # Flip order for better visualization in pymol
        # denoised
        denoised_xyz_stack = torch.stack(denoised_xyz_stack)
        self.denoised_xyz_stack = torch.flip(denoised_xyz_stack, [0,],)
        # px0 
        px0_xyz_stack = torch.stack(px0_xyz_stack)
        self.px0_xyz_stack = torch.flip(px0_xyz_stack, [0,],)
        # For logging -- don't flip
        self.plddt_stack = torch.stack(plddt_stack)

        # final seq
        self.final_seq = seq_stack[-1]
        # Output glycines, except for motif region
        # padding by 7, representing glycine
        self.final_seq = torch.where(
            torch.argmax(seq_init, dim=-1) == 21, 7, torch.argmax(seq_init, dim=-1)
        )  
        aa_seq = seq2chars(self.final_seq)
        print(aa_seq)

        # bfacts: from final_seq
        self.bfacts = torch.ones_like(self.final_seq.squeeze())
        # print(self.bfacts)
        # make bfact=0 for diffused coordinates
        self.bfacts[torch.where(torch.argmax(seq_init, dim=-1) == 21, True, False)] = 0

    def to_trb(self, out_prefix:str, start_time):
        '''
        need self.plddt_stack
        create *.trb
        '''
        trb = dict(
            config=OmegaConf.to_container(self.sampler._conf, resolve=True),
            plddt=self.plddt_stack.cpu().numpy(),
            device=self.device,
            time = time.time() - start_time,
        )
        if hasattr(self.sampler, "contig_map"):
            for key, value in self.sampler.contig_map.get_mappings().items():
                trb[key] = value
        # export
        outfile = f"{out_prefix}.trb"
        with open(outfile, "wb") as f:
            pickle.dump(trb, f)
        print("meta data:", outfile)

    def from_trb(self, num:int):
        trb_file = f"{self.sampler._conf.inference.output_prefix}_{num}.trb"
        with open(trb_file, 'rb') as f:
            return pickle.load(f)

    def display_pdb(self, num:int, cartoon=None):
        '''
        display 3D given *.pdb
        '''
        # pdb file
        pdb_file = f"{self.sampler._conf.inference.output_prefix}_{num}.pdb"
        print(pdb_file)
        pdb_str = open(pdb_file,'r').read()
        # parameters
        params = {
            'hbondCutoff': 4.0
        }
        default_cartoon = {
            'colorscheme': {
                'prop':'b',
                'gradient':'roygb',
                'min':0,
                'max':100
            },
            'color': 'spectrum',
        }
        cartoon = cartoon if cartoon else default_cartoon

        # display
        view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js')
        view.addModel(pdb_str, 'pdb', params)
        view.setStyle({"model":0}, {'cartoon':cartoon})
        view.zoomTo()
        view.show()
        return pdb_file