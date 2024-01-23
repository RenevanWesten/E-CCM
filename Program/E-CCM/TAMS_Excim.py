#TAMS on the E-CCM

import os
import numpy as np
from tqdm import trange
import argparse

from TAMS_full import TAMS
from MOC_model import MOCModel

# We are going to run the TAMS algorithm on the MOC model
# Range of parameters: E_A in [0, 0.5] Sv, noise in [0, 0.5]
# 100 points in each direction
# TAMS is run 15 times for each parameter set, with 1000 ensemble members, for 100 model years
# At each iteration of TAMS, we discard nc=10 ensemble members (1% of all ensemble members)

def init_TAMS(tams, C, start, z=0):
    if start == "on":
        tams.set_check_on(C.is_on)
        tams.set_check_off(C.check_F_collapse)
        tams.set_traj_func(C.run, stop_off=True)
        init = C.on
    else:
        tams.set_check_on(C.check_S_collapse)
        tams.set_check_off(C.is_on)
        tams.set_traj_func(C.run, stop_on=True)
        init = C.off
    tams.set_score(C.score, z=z, on=(start=="on"))
    return tams, init

parser = argparse.ArgumentParser(description='Run TAMS on the MOC model')
parser.add_argument('--N_Ea', type=int, default=100, help='Number of points in the E_A direction')
parser.add_argument('--N_noise', type=int, default=100, help='Number of points in the noise direction')
parser.add_argument('--Ea', type=float, help='Specific value of E_A')
parser.add_argument('--noise', type=float, help='Specific value of noise')
parser.add_argument('--nc', type=int, default=10, help='Number of ensemble members discarded at each iteration')
parser.add_argument('--nb_runs', type=int, default=15, help='Number of runs for each parameter set')
parser.add_argument('--N', type=int, default=1000, help='Number of ensemble members')
parser.add_argument('--tmax', type=int, default=100, help='Number of model years')
parser.add_argument('--start', type=str, default='on', help='Starting state of the model, type "on" or "off"')
parser.add_argument('--seed', type=int, default=0, help='Seed for the model, default is 0')
parser.add_argument('--delta_t', type=int, default=20*86400, help='Time step of the model, default is 10 days')
parser.add_argument('--temp_params', type=float, help="Value of lambda")
parser.add_argument('--ice_params', type=float, nargs='+', help="Values of f_ice q_factor")
parser.add_argument('--output_dir', type=str, default='./', help='Name of the output directory')
parser.add_argument('--min_proba', type=float, default=1e-9, help='Minimum probability to be computed by TAMS')
parser.add_argument('--min_level', type=float, default=0, help='Minimum target level of the score function')
args = parser.parse_args()

N_Ea, N_noise = args.N_Ea, args.N_noise
nc = args.nc
nb_runs, N = args.nb_runs, args.N
tmax = args.tmax
start = args.start

if args.output_dir[-1] != "/":
    args.output_dir += "/"

tams = TAMS(N, nc, tmax)


if args.temp_params is not None:
    if args.ice_params is not None:
        C = MOCModel(delta_t=args.delta_t, seed=args.seed, temp_params=args.temp_params, ice_params=args.ice_params)
    else:
        C = MOCModel(delta_t=args.delta_t, seed=args.seed, temp_params=args.temp_params)
else:
    C = MOCModel(delta_t=args.delta_t, seed=args.seed)

#Get the limits of the on and off state
E_min_off_branch	= 240 - C.E_A_index_off
C.set_E_A(0.5)
E_max_on_branch		= C.E_A_index_on

#Now convert to the original data
E_A 			    = np.linspace(0.0, 0.5, 101)
E_min_off_branch	= E_A[E_min_off_branch]
E_max_on_branch		= E_A[E_max_on_branch]

if args.Ea is not None and args.noise is not None:

    if args.temp_params is not None:
        if args.ice_params is not None:
            name = "TAMS_temp_ice_Ea_{}_noise_{}_start_{}_nc_{}_tmax_{}".format(C.str_param(args.Ea),
                                                                                C.str_param(args.noise),
                                                                                start, nc, tmax)
        else:
            name = "TAMS_temp_Ea_{}_noise_{}_start_{}_nc_{}_tmax_{}".format(C.str_param(args.Ea),
                                                                             C.str_param(args.noise),
                                                                             start, nc, tmax)
    else:
        name = "TAMS_Ea_{}_noise_{}_start_{}_nc_{}_tmax_{}".format(C.str_param(args.Ea),
                                                                   C.str_param(args.noise),
                                                                   start, nc, tmax)

    C.set_E_A(args.Ea)
    C.set_noise(args.noise)

    tams, init = init_TAMS(tams, C, start)

    p = np.zeros(nb_runs)
    for i in trange(nb_runs):
        tams.reset_seed(i)
        p[i], _, _, traj = tams.run(init)
    print("Mean transition probability: {}".format(np.mean(p)))

    np.savez(args.output_dir + name, probas=p, trajectories=traj)

else:

    if args.temp_params is not None:
        if args.ice_params is not None:
            name = "TAMS_temp_ice_start_{}_nc_{}_tmax_{}".format(start, nc, tmax)
        else:
            name = "TAMS_temp_start_{}_nc_{}_tmax_{}".format(start, nc, tmax)
    else:
        name = "TAMS_start_{}_nc_{}_tmax_{}".format(start, nc, tmax)

    tams.set_kmax(args.min_proba)

    E_A   = np.linspace(0.0, 0.5, N_Ea)
    noise = np.linspace(0, 0.5, N_noise)

    if name + ".npz" not in os.listdir(args.output_dir):
        np.savez(args.output_dir + name, probas=np.full((N_noise,N_Ea,nb_runs), np.nan))
    
    with np.load(args.output_dir + name + ".npz") as data:
        all_probas = data['probas']

    for i, ea in enumerate(E_A[::-1]):
        for j, n in enumerate(noise[::-1]):

            print("E_A = {}, noise = {}".format(ea, n))

            if ea < E_min_off_branch or ea > E_max_on_branch:
                #Outside the steady state regimes, transition probabilities are meaningless
                all_probas[N_noise-j-1, N_Ea-i-1, :] = np.zeros(nb_runs) -1
                np.savez(args.output_dir + name, probas=all_probas)
                continue

            if np.any(np.isnan(all_probas[N_noise-j-1, N_Ea-i-1,:])) or np.any(all_probas[N_noise-j-1, N_Ea-i-1,:] < args.min_proba):  

                C.set_E_A(ea)
                C.set_noise(n)
                C.reset_seed(args.seed)
                tams, init = init_TAMS(tams, C, start, z=args.min_level)

                #Get the probability from the general array, usually nan
                p = all_probas[N_noise-j-1, N_Ea-i-1,:]

                for k in trange(nb_runs):
                    #Loop over all the TAMS iterations
                    print(p)

                    if np.isnan(p[k]) == False and p[k] > args.min_proba:
                        #Probability is already determined, continue to next TAMS iteration
                        continue

                    #Conduct TAMS
                    tams.reset_seed(k)
                    p[k], _, _, _ = tams.run(init)

                    if p[k] < args.min_proba:
                        #Probability is too low, set to very small value
                        p[k]	= 10**(-30.) 

                    #Save TAMS after every succesfull determination
                    all_probas[N_noise-j-1, N_Ea-i-1, :] = p
                    np.savez(args.output_dir + name, probas=all_probas)

                    if p[k] < args.min_proba:
                        #Probability is too low, set to very small value
                        break
