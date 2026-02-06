import numpy as np

import vmc

def e_l_alpha(R:np.ndarray, alpha: float):
  R_norm = np.linalg.vector_norm(R)
  Rn_invs = 1/R_norm
  return -Rn_invs - 0.5*alpha*(alpha - 2*Rn_invs)

def psi_alpha(R:np.ndarray, alpha:float):
  R_norm = np.linalg.vector_norm(R)
  return np.exp(-alpha * R_norm)

def logder(R:np.ndarray, alpha:float):
  R_norm = np.linalg.vector_norm(R)
  return -R_norm

if __name__=="__main__":
  vmc.variational_mc(1, 3, psi_alpha, e_l_alpha, logder, "hydrogen atom")