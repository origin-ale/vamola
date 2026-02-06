import numpy as np
import time

import config_walk as cw
import vmc

import harmonic_oscillator
import hydrogen
import helium

def test_vmc(particles: int, 
                   dims: int, 
                   psi_alpha: function, 
                   e_l_alpha: function, 
                   alphas: list,
                   syst_name: str = 'unknown system'):
  """
  Returns energies and standard deviations obtained by variational Monte Carlo algorithm for a list of alphas.\\
  All function parameters should be over config space and depend on one parameter alpha.
  
  Parameters
  --
  particles: int
    The number of particles in the system
  dims: int
    The number of spatial dimensions
  psi_alpha: function(R, alpha)
    The variational approximation for the wave function.
  e_l_alpha: function(R, alpha)
    The local energy of the system. Implements the hamiltonian and depends on the variational WF.\\
    Calculated analytically.
  alphas: list[float]
    The values of the variational parameter to calculate the energy for.
  """

if __name__ == "__main__":
  st = time.perf_counter()
  alphas_ho = [0.4, 0.45, 0.5, 0.55, 0.6]
  for a in alphas_ho:
    energy, _, le_stdev, _ = vmc.vmc_iteration(0,
                                                             a, 
                                                             harmonic_oscillator.e_l_alpha, 
                                                             harmonic_oscillator.psi_alpha, 
                                                             harmonic_oscillator.logder, 
                                                             1, 
                                                             1,
                                                             40, 
                                                             3000, 
                                                             400,
                                                             1, 
                                                             st)
    print(f"alpha = {a}\tE = {energy}\tvar={le_stdev}")