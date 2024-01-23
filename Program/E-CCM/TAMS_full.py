#TAMS procedure

import numpy as np
import matplotlib.pyplot as plt

class TAMS():
    
    def __init__(self, N, nc, Tmax, seed=None):
        self.N, self.nc, self.Tmax = N, nc, Tmax
        self.kmax = np.inf
        self.rng = np.random.default_rng(seed=seed)
        
    def reset_seed(self, seed):
        self.rng.bit_generator.state = np.random.PCG64(seed).state
    
    def set_kmax(self, q_min):
        """
        Set the maximum number of iterations given a minimum probability 
        """
        self.kmax = np.ceil(np.log(q_min)/np.log(1-self.nc/self.N))

    def set_score(self, score_function, *fixed_args, **fixed_kwargs):
        self.comp_score = lambda traj : score_function(traj, *fixed_args, **fixed_kwargs)
    
    def set_traj_func(self, traj_function, *fixed_args, **fixed_kwargs):
        if 'need_idx' in fixed_kwargs:
            fixed_kwargs = {key: value for key, value in fixed_kwargs.items() if key != 'need_idx'}
            self.comp_traj = lambda n, t, ic, idx : traj_function(n, t, ic, idx, *fixed_args, **fixed_kwargs)
        else:
            self.comp_traj = lambda n, t, ic : traj_function(n, t, ic, *fixed_args, **fixed_kwargs)
            
    def set_check_on(self, on_function, *fixed_args, **fixed_kwargs):
        self.check_on = lambda traj : on_function(traj, *fixed_args, **fixed_kwargs)
        
    def set_check_off(self, off_function, *fixed_args, **fixed_kwargs):
        self.check_off = lambda traj : off_function(traj, *fixed_args, **fixed_kwargs)
    
    def run(self, ic, zmin=0, zmax=1, need_idx=False):

        k, w = 0, 1
        
        if not need_idx:          
            traj = self.comp_traj(self.N, self.Tmax, ic)
        else:
            traj = self.comp_traj(self.N, np.zeros(self.N,dtype=int), ic, np.arange(self.N))
            all_ind = np.arange(self.N)
            
        Nt = traj.shape[1]
        nb_total_timesteps = Nt*self.N
        
        score = self.comp_score(traj)
        onzone = self.check_on(traj)
        offzone = self.check_off(traj)

        score[onzone], score[offzone] = zmin, zmax

        first_idx_off = np.argmax(offzone, axis=1)
        for i in range(self.N):
            if first_idx_off[i]>0:
                score[i,first_idx_off[i]+1:] = np.nan
                traj[i,first_idx_off[i]+1:] = np.nan

        Q = np.nanmax(score,axis=1)

        while len(np.unique(Q)) > self.nc and k<self.kmax:
            
            threshold = np.unique(Q)[self.nc-1] #Because Python counts from 0
            idx, other_idx = np.nonzero(Q<=threshold)[0], np.nonzero(Q>threshold)[0]
            Q_min = Q[idx]
            #print(Q_min)

            #Update weights
            w *= (1-len(idx)/self.N)
            
            # Create new trajectories
            new_ind = self.rng.choice(other_idx, size=len(idx))
            restart = np.argmax(score[new_ind]>=Q_min[:,np.newaxis], axis=1)
            init_clone = traj[new_ind,restart]
            
            all_length = Nt - restart
            nb_total_timesteps += np.sum(all_length)

            if not need_idx:            
                new_traj = self.comp_traj(len(idx), self.Tmax*np.max(all_length)/Nt, init_clone)
            else:
                new_traj = self.comp_traj(len(idx), restart, init_clone, all_ind[new_ind])
                all_ind[idx] = all_ind[new_ind]
            
            new_score = self.comp_score(new_traj)
            
            for i in range(len(idx)):
                t, r, l = idx[i], restart[i], all_length[i]
                
                traj[t,:r] = traj[new_ind[i],:r]
                traj[t,r:] = new_traj[i,:l]
                
                score[t,:r+1] = score[new_ind[i],:r+1]
                score[t,r+1:] = new_score[i,1:l]
                
                onzone = self.check_on(traj[t])
                offzone = self.check_off(traj[t])
                score[t,onzone], score[t,offzone] = zmin, zmax
                
                first_idx_off = np.argmax(offzone)
                if first_idx_off>0:
                    score[t,first_idx_off+1:] = np.nan
                    traj[t,first_idx_off+1:] = np.nan
                
            k += 1
            Q = np.nanmax(score,axis=1)

        if not need_idx:
            return w*np.count_nonzero(Q>=1)/self.N, k, nb_total_timesteps, traj
        return w*np.count_nonzero(Q>=1)/self.N, k, nb_total_timesteps, all_ind, traj

#%%
