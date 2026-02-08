import time
import numpy as np
from functools import partial

import config_walk as cw
import parser

def update_alpha(alpha: float, 
                 energy: float, 
                 logder:function, 
                 e_l_alpha:function, 
                 samples: np.ndarray, 
                 gamma = 1.):
  """"
  Use damped steepest descent to calculate a new value of the variational parameter alpha,\\
  evaluating dE/d(alpha) from VMC samples and local energy (Thijssen eq. 12.13).

  Parameters
  ---
  alpha: float
    The current value of alpha
  energy: float
    The current value of the energy
  logder: function(R, alpha) -> float
    d/d(alpha) [ln psi_alpha], used in varying alpha.\\
    Calculated analytically.
  e_l_alpha: function(R, alpha) -> float
    The local energy of the system. Implements the hamiltonian and depends on the variational WF.\\
    Calculated analytically.
  samples: array of configurations
    The samples to average over
  gamma: float
    The steepest descent parameter
  
  Returns
  ---
  The new value of alpha
  """
  t1  = cw.sample_avg(samples, lambda x: e_l_alpha(x,alpha)*logder(x, alpha))
  t2 = -energy * cw.sample_avg(samples, lambda x:logder(x, alpha))
  der = 2 * (t1 + t2)
  return alpha - gamma * der

def vmc_sample(alpha:float,          
               psi_alpha: function,
               particles: int,
               dims: int, 
               walker_step: float, 
               walker_n: int, 
               steps: int, 
               thermal: int, 
               print_interval:int,
               st:float):
  psi = partial(psi_alpha, alpha=alpha)
  walkers = [cw.ConfigWalker(particles,dims,walker_step) for i in range(0,walker_n)]
  samples = []
  elapsed_time = time.perf_counter() - st
  next_print = int(elapsed_time) + print_interval
  for i in range(0, steps):
    for w in walkers:
      w.metropolis_step(psi)
      if i >= thermal:
        samples.append(w.current_config())
    elapsed_time = time.perf_counter() - st
    if int(elapsed_time) % print_interval == 0 and int(elapsed_time) >= next_print: 
      print(f"{elapsed_time:.1f} s: Metropolis step {i} ({i/steps * 100:.1f}%)")
      next_print = int(elapsed_time) + print_interval
    
  return samples
  
def vmc_energy(alpha: float, e_l_alpha: np.ndarray, samples: list, st:float):
  e_l = partial(e_l_alpha, alpha=alpha)
  elapsed_time = time.perf_counter() - st
  print(f"{elapsed_time:.1f} s: calculating energy and stdev from {len(samples)} samples...")
  energy = cw.sample_avg(samples, e_l)
  le_stdev = cw.sample_stdev(samples, e_l, energy)
  return energy, le_stdev


def variational_mc(particles: int, 
                   dims: int, 
                   walker_step: float,
                   psi_alpha: function, 
                   e_l_alpha: function, 
                   logder: function, 
                   syst_name: str = 'unknown system'):
  """
  Run a variational Monte Carlo algorithm.\\
  All function-type inputs should be over config space and depend on one parameter alpha.
  
  Parameters
  --
  particles: int
    The number of particles in the system
  dims: int
    The number of spatial dimensions
  psi_alpha: function(R, alpha) -> float
    The variational approximation for the wave function.
  e_l_alpha: function(R, alpha) -> float
    The local energy of the system. Implements the hamiltonian and depends on the variational WF.\\
    Calculated analytically.
  logder: function(R, alpha) -> float
    d/d(alpha) [ln psi_alpha], used in varying alpha.\\
    Calculated analytically.

  Returns
  ---
  alphas: list of float
    The values of alpha considered
  energies: list of float
    The energies for each alpha considered until convergence
  stdevs: list of float
    The energy stdevs for each alpha
  """
  cl_args = parser.vmc_args()
  start_alpha = cl_args.alpha
  walker_n = cl_args.walkers
  steps = cl_args.steps
  thermal = cl_args.thermalization
  print_interval = cl_args.print
  threshold = cl_args.convergence

  alphas = []
  energies = []
  stdevs = []

  energy = 0.
  le_stdev = np.inf
  energy_diff = np.inf
  alpha = start_alpha
  st = time.perf_counter()
  j = 1
  print("Running VMC for " + f"{syst_name}.")
  print(f"Walkers: {walker_n}")
  print(f"Metropolis steps: {steps}\t\tThermalization: {thermal}\t\tInfo print interval: {print_interval} s")
  print(f"Starting alpha: {alpha}\t\tEnergy convergence threshold: {threshold}")
  while energy_diff > threshold:
    print(8*'-' + f" Iteration {j} " + 32*'-')
    samples = vmc_sample(alpha, psi_alpha, particles, dims, walker_step, walker_n, steps, thermal, print_interval, st)

    energy_last = energy
    energy, le_stdev = vmc_energy(alpha, e_l_alpha, samples, st)
    elapsed_time = time.perf_counter() - st
    energy_diff = abs(energy_last-energy)
    
    alphas.append(alpha)
    energies.append(energy)
    stdevs.append(le_stdev)

    if energy_diff > threshold:
      print(f"{elapsed_time:.1f} s: calculated energy = {energy} with stdev = {le_stdev}. Now updating alpha...")
      alpha = update_alpha(alpha, energy, logder, e_l_alpha, samples)
      elapsed_time = time.perf_counter() - st
      print(f"{elapsed_time:.1f} s: set alpha = {alpha}")

    j+=1

  print(8*'=' + f" Converged with alpha = {alpha} " + 64*'=')
  print(f"Energy = {energy} with stdev = {le_stdev}")

  return alphas, energies, stdevs

  