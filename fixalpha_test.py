import numpy as np
import matplotlib.pyplot as plt
import time

import vmc

import harmonic_oscillator as ho
import hydrogen as hy
import helium as he

def wf_energy(alpha, psi_alpha, e_l_alpha):
  x = np.linspace(0, 5, 1000)

  wf_weights = psi_alpha(x, alpha) ** 2
  e_l_x = e_l_alpha(x, alpha)
  energy = np.average(e_l_x, weights= wf_weights)
  return energy

def test_vmc(particles: int, 
            dims: int, 
            walker_step: int,
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
  print("Running VMC energy test for " + syst_name + "...")
  energies = []
  for a in alphas:
    st = time.perf_counter()
    samples = vmc.vmc_sample(a, 
                            psi_alpha, 
                            particles, 
                            dims,
                            walker_step,
                            10, 
                            100000, 
                            400,
                            1, 
                            st)
    r1_samples = [np.linalg.norm(s[0:2]) for s in samples]
    sample_hist = plt.hist(r1_samples, bins = 1000)
    plt.show()
    energy, le_stdev = vmc.vmc_energy(a, e_l_alpha, samples)
    energies.append(energy)
    print(f"alpha = {a}\tE = {energy}\tvar={le_stdev}")
  return energies

if __name__ == "__main__":
  # ho_step = .33
  # alphas_ho = [0.4, 0.45, 0.5, 0.55, 0.6]
  # energies_ho = test_vmc(1,1,ho_step, ho.psi_alpha, ho.e_l_alpha, alphas_ho, "harmonic oscillator")

  # hy_step = 1.5
  # alphas_hy = [0.8, 0.9, 1, 1.1, 1.2]
  # energies_hy = test_vmc(1,3,hy_step,hy.psi_alpha, hy.e_l_alpha, alphas_hy, "hydrogen atom")

  he_step = 2
  alphas_he = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]
  energies_he = test_vmc(2,3,he_step,he.psi_alpha, he.e_l_alpha, alphas_he, "helium atom")

  fig, axs = plt.subplots(1,3, figsize = [19.2, 3.6])

  # axs[0].plot(alphas_ho, energies_ho)
  # f = lambda x: 0.5*x + 1/(8*x)
  # energies_expected=[f(a) for a in alphas_ho]
  # axs[0].plot(alphas_ho, energies_expected)

  # axs[1].plot(alphas_hy, energies_hy)

  axs[2].plot(alphas_he, energies_he)
  
  plt.show()