import numpy as np

import vmc
import output_utils as ou

def e_l_alpha(x:float, alpha: float):
  return alpha + x**2 * (0.5 - 2*(alpha**2))

def psi_alpha(x:float, alpha:float):
  return np.exp(-alpha * (x**2))

def logder(x:float, alpha:float):
  return -(x**2)

if __name__=="__main__":
  walker_step = .4
  alphas, energies, stdevs = vmc.variational_mc(1, 
                                                1, 
                                                walker_step, 
                                                psi_alpha, 
                                                e_l_alpha, 
                                                logder, 
                                                "1D harmonic oscillator")
  iterations = [i+1 for i in range(0, len(alphas))]
  ou.lists_to_file(ou.conv_name("ho"), iterations, alphas, energies, stdevs, headers=["iteration", "alpha", "energy", "stdev"])