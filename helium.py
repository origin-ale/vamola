import numpy as np

import vmc
import output_utils as ou

def e_l_alpha(R:np.ndarray, alpha: float):
  r1 = R[0:3]
  r2 = R[3:6]
  r12 = r1-r2
  r1_norm = np.linalg.vector_norm(r1)
  r2_norm = np.linalg.vector_norm(r2)
  rv1  = r1/r1_norm
  rv2 = r2/r2_norm
  r12_norm = np.linalg.vector_norm(r12)
  if r12_norm == 0: return 0

  x = 1 + alpha * r12_norm

  o0 = -4
  o1 = 1/r12_norm
  o2 = np.dot(rv1 - rv2,r1-r2) /(r12_norm * x**2)
  o3 = - 1/(r12_norm * x**3)
  o4 = - 1/(4 * x**4)

  return o0 + o1 + o2 + o3 + o4

def psi_alpha(R:np.ndarray, alpha:float):
  r1 = R[0:3]
  r2 = R[3:6]
  r12 = r1-r2
  r1_norm = np.linalg.vector_norm(r1)
  r2_norm = np.linalg.vector_norm(r2)
  r12_norm = np.linalg.vector_norm(r12)

  return np.exp(-2*r1_norm - 2*r2_norm) * np.exp(r12_norm/(2*(1+alpha*r12_norm)))
  

def logder(R:np.ndarray, alpha:float):
  r1 = R[0:3]
  r2 = R[3:6]
  r12 = r1-r2
  r12_norm = np.linalg.vector_norm(r12)
  return - 0.5 * (r12_norm/(1 + alpha * r12_norm))**2


if __name__=="__main__":
  walker_step = .4
  alphas, energies, stdevs = vmc.variational_mc(2, 
                                                3, 
                                                walker_step, 
                                                psi_alpha, 
                                                e_l_alpha, 
                                                logder, 
                                                "helium atom")
  iterations = [i+1 for i in range(0, len(alphas))]
  ou.lists_to_file(ou.conv_name("he"), iterations, alphas, energies, stdevs, headers=["iteration", "alpha", "energy", "stdev"])