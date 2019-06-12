import os, sys

import numpy as np

from tqdm import tqdm_notebook as tqdm

# import solver related modules

# add path to MCEq installation
#sys.path.append("/home/eganster/software/MCEq")

# import solver related modules
from MCEq.core import MCEqRun
from mceq_config import config


# update some paths in the mceq_config

# Directory where the data files for the calculation are stored
#config["data_dir"] = os.path.join("/home/eganster/software/MCEq", "data")
# full path to libmkl_rt.[so/dylib] (only if kernel=='MKL')
config["MKL_path"] = os.path.join("/home/eganster/.local/lib", "libmkl_rt.so")

# define function to access MCEq solutions more easily

def get_solutions(mceq_run, particle_ids, mag=3):
        return_fluxes = {}
        for pid in particle_ids:
            return_fluxes[pid] = mceq_run.get_solution(pid, mag)

        return return_fluxes
    

# set debug level for MCEq (default = 1)
config["debug_level"] = 0


def run_MCEq(primary_model,
             interaction_model="SIBYLL2.3c",
             density_profiles=[("MSIS00_IC", ("SouthPole", "January")),
                                ("MSIS00_IC", ("SouthPole", "July")),
                              ],
             particle_ids=["total_numu", "total_antinumu"],
             cosz_lim=[-1.0, 1.0],
             cosz_steps=50,
             emag=3,
            ):

    # define equidistant grid in cos(theta) with cosz_steps steps
    theta_grid = np.arccos(np.linspace(cosz_lim[0], cosz_lim[1], cosz_steps))
    theta_grid *= 180.0/np.pi
    
    # temporarily result dict
    flux_for_density = {}
    
    ## loop over all density profiles
    #for density in density_profiles:
    for density in tqdm(density_profiles, desc="Density"):
        print "="*60
        print "Current atmosphere model:", density[0], "--", density[1][0], density[1][1]
        print "-"*60
        
        # set atmosphere model string for result dictionary
        if density[1][1] is not None:
            density_str = density[0]+density[1][0]+density[1][1]
        else:
            density_str =  density[0]+density[1][0]
        
        # update mceq_config with the current atmosphere model
        config["density_model"] = density
        
        # create instance of MCEqRun class
        mceq_run = MCEqRun(interaction_model = interaction_model,
                           primary_model     = primary_model,
                           theta_deg         = 0.0,        # updated later
                           **config)
        
        # obtain energy grid (fixed) of the solution for the x-axis of the plots
        e_grid = mceq_run.e_grid
        
        # update dictionary
        flux_for_density[density_str] = {}
        for flux_str in particle_ids:
            flux_for_density[density_str][flux_str] = np.zeros((len(theta_grid), len(e_grid)))

        ## loop over all theta bins
        #for theta_id, theta in enumerate(theta_grid):
        for theta_id, theta in enumerate(tqdm(theta_grid, desc="Theta")):
            print "-"*60
            print "Current theta:", theta
            
            # Set/update the zenith angle
            mceq_run.set_theta_deg(theta)
            # Run the solver
            mceq_run.solve()
            
            # get fluxes
            flux_solutions = get_solutions(mceq_run, particle_ids, mag=emag)
            
            # store fluxes in result dictionary'
            for flux_str in particle_ids:
                flux_for_density[density_str][flux_str][theta_id,:] = flux_solutions[flux_str]
                #print flux_for_density[density_str][flux_str][theta_id,:]
        
    # average density models:
    fluxes = {}
    for flux_str in particle_ids:
        fluxes[flux_str] = np.zeros((len(theta_grid), len(e_grid)))

    for flux_str in particle_ids:
        for density in flux_for_density.keys(): # loop over all density models
            fluxes[flux_str] += flux_for_density[density][flux_str]
        fluxes[flux_str] /= len(flux_for_density.keys())*1.0
        
    # add e_grid and theta_grid
    #fluxes["e_grid"] = e_grid
    #fluxes["theta_grid"] = theta_grid
    
    return fluxes