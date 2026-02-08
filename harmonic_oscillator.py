import numpy as np

import vmc
import output_utils as ou

def e_l_alpha(R:np.ndarray, alpha: float):
  x = np.linalg.vector_norm(R)
  return alpha + x**2 * (0.5 - 2*(alpha**2))

def psi_alpha(R:np.ndarray, alpha:float):
  x = np.linalg.vector_norm(R)
  return np.exp(-alpha * (x**2))

def logder(R:np.ndarray, alpha:float):
  x = np.linalg.vector_norm(R)
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