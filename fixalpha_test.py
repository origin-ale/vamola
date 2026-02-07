import numpy as np
import matplotlib.pyplot as plt
import time

import config_walk as cw
import vmc

import harmonic_oscillator as ho
import hydrogen as hy
import helium as he

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
  print("Running test for " + syst_name + "...")
  energies = []
  for a in alphas:
    st = time.perf_counter()
    samples = vmc.vmc_sample(a, 
                            psi_alpha, 
                            particles, 
                            dims,
                            10, 
                            10000, 
                            5000,
                            1, 
                            st)
    energy, le_stdev = vmc.vmc_energy(a, e_l_alpha, samples)
    energies.append(energy)
    # wf_weights = psi_alpha(x, a) ** 2
    # avg_le = np.average(e_l,weights=wf_weights)
    # avg_les.append(avg_le)
    print(f"alpha = {a}\tE = {energy}\tvar={le_stdev}")
  return energies

if __name__ == "__main__":
  alphas_ho = [0.4, 0.45, 0.5, 0.55, 0.6]
  energies_ho = test_vmc(1,1,ho.psi_alpha, ho.e_l_alpha, alphas_ho, "harmonic oscillator")

  alphas_hy = [0.8, 0.9, 1, 1.1, 1.2]
  energies_hy = test_vmc(1,3,hy.psi_alpha, hy.e_l_alpha, alphas_hy, "hydrogen atom")

  alphas_he = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25]
  energies_he = test_vmc(2,3,he.psi_alpha, he.e_l_alpha, alphas_he, "helium atom")

  fig, axs = plt.subplots(1,3, figsize = [19.2, 3.6])
  axs[0].plot(alphas_ho, energies_ho)
  f = lambda x: 0.5*x + 1/(8*x)
  energies_expected=[f(a) for a in alphas_ho]
  axs[0].plot(alphas_ho, energies_expected)

  axs[1].plot(alphas_hy, energies_hy)
  # f = lambda x: 0.5*x + 1/(8*x)
  # energies_expected=[f(a) for a in alphas_hy]
  # axs[1].plot(alphas_hy, energies_expected)

  axs[2].plot(alphas_he, energies_he)
  
  plt.show()