import numpy as np

import config_walk as cw
import vmc

def e_l_alpha(x:float, alpha: float):
  return alpha + x**2 * (0.5 - 2*alpha**2)

def psi_alpha(x:float, alpha:float):
  return np.exp(-alpha * x**2)

def logder(x:float):
  return -x**2

if __name__=="__main__":
  vmc.variational_mc(psi_alpha, e_l_alpha, logder, "1D harmonic oscillator")