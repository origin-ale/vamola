import numpy as np
import matplotlib.pyplot as plt
import time

import parser
import vmc
import output_utils as ou

import harmonic_oscillator as ho
import hydrogen as hy
import helium as he

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
  alphas: list of float
    The values of the variational parameter to calculate the energy for.
  """
  print("Running VMC energy test for " + syst_name + "...")
  energies = []
  stdevs = []
  st = time.perf_counter()
  for a in alphas:
    samples = vmc.vmc_sample(a, 
                            psi_alpha, 
                            particles, 
                            dims,
                            walker_step,
                            400, 
                            30000, 
                            4000,
                            5, 
                            st)
    energy, le_stdev = vmc.vmc_energy(a, e_l_alpha, samples, st)
    energies.append(energy)
    stdevs.append(le_stdev)
    elapsed_time = time.perf_counter() - st
    print(f"{elapsed_time:.1f} s: calculated alpha = {a}\tE = {energy}\tvar={le_stdev}")
  return energies,stdevs

if __name__ == "__main__":
  cli_args = parser.test_args()
  systs = cli_args.systems.lower().split()
  syst_n = 0

  if "ho" in systs:
    ho_step = .33
    alphas_ho = [0.4, 0.45, 0.5, 0.55, 0.6]
    energies_ho, stdevs_ho = test_vmc(1,1,ho_step, ho.psi_alpha, ho.e_l_alpha, alphas_ho, "harmonic oscillator")
    ou.lists_to_file(ou.fa_name("ho"), alphas_ho, energies_ho, stdevs_ho, headers=["alpha", "energy", "stdev"])
    syst_n += 1

  if "hy" in systs:
    hy_step = 1.5
    alphas_hy = [0.8, 0.9, 1, 1.1, 1.2]
    energies_hy, stdevs_hy = test_vmc(1,3,hy_step,hy.psi_alpha, hy.e_l_alpha, alphas_hy, "hydrogen atom")
    ou.lists_to_file(ou.fa_name("hy"), alphas_hy, energies_hy, stdevs_hy, headers=["alpha", "energy", "stdev"])
    syst_n += 1

  if "he" in systs:
    he_step = 2
    alphas_he = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225]
    energies_he, stdevs_he = test_vmc(2,3,he_step,he.psi_alpha, he.e_l_alpha, alphas_he, "helium atom")
    ou.lists_to_file(ou.fa_name("he"), alphas_he, energies_he, stdevs_he, headers=["alpha", "energy", "stdev"])
    syst_n += 1


  if syst_n > 1: fig, axs = plt.subplots(1,syst_n, figsize = [6.4*syst_n, 3.6])

  if syst_n == 1 :
    fig,ax = plt.subplots(1,1, figsize = [6.4, 3.6])
    axs = [ax,]

  current_ax = 0
  if "ho" in systs:
    axs[current_ax].plot(alphas_ho, energies_ho, label = "VMC, harm. osc.")
    exp_ho = [0.51241, 0.502764, 0.5, 0.502326, 0.50841]
    axs[current_ax].plot(alphas_ho, exp_ho, label = "expected, harm. osc.")
    axs[current_ax].legend()
    current_ax+=1


  if "hy" in systs: 
    axs[current_ax].plot(alphas_hy, energies_hy, label = "VMC, H")
    exp_hy = [-0.47962, -0.49491, -0.5, -0.49512, -0.48013]
    axs[current_ax].plot(alphas_hy, exp_hy, label = "expected, H")
    axs[current_ax].legend()
    current_ax+=1

  if "he" in systs: 
    axs[current_ax].plot(alphas_he, energies_he, label = "VMC, He")
    exp_he = [-2.87134, -2.87534, -2.87703, -2.87804, -2.87783, -2.87813, -2.87674, -2.87461]
    axs[current_ax].plot(alphas_he, exp_he, label = "expected, He")
    axs[current_ax].legend()
    current_ax+=1
  
  plt.show()