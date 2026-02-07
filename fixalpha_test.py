import numpy as np
import matplotlib.pyplot as plt
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
  # alphas_ho = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
  energies_ho = []
  x = np.linspace(-10,10, 20000)
  e_ls=[]
  avg_les = []
  for a in alphas_ho:
    energy, _, le_stdev, _ = vmc.vmc_iteration(0,
                                              a, 
                                              harmonic_oscillator.e_l_alpha, 
                                              harmonic_oscillator.psi_alpha, 
                                              harmonic_oscillator.logder, 
                                              1, 
                                              1,
                                              10, 
                                              10000, 
                                              5000,
                                              1, 
                                              st)
    e_l = harmonic_oscillator.e_l_alpha(x, a)
    energies_ho.append(energy)
    wf_weights = harmonic_oscillator.psi_alpha(x, a) ** 2
    avg_le = np.average(e_l,weights=wf_weights)
    avg_les.append(avg_le)
    e_ls.append(e_l)
    print(f"alpha = {a}\tE = {energy}\tvar={le_stdev}; average of local energy {avg_le}")
  

  plt.plot(alphas_ho, energies_ho)
  f = lambda x: 0.5*x + 1/(8*x)
  energies_expected=[f(a) for a in alphas_ho]
  plt.plot(alphas_ho, energies_expected)
  
  plt.show()