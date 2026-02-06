import time
import numpy as np
from functools import partial

import config_walk as cw
import parser

def update_alpha(alpha: float, energy: float, logder:function, e_l_alpha:function, samples: np.ndarray, gamma = 1.):
  der = 2 * (cw.sample_avg(samples, lambda x: e_l_alpha(x,alpha)*logder(x, alpha)) - energy * cw.sample_avg(samples, lambda x:logder(x, alpha)))
  return alpha - gamma * der

def vmc_iteration(energy, alpha, e_l_alpha, psi_alpha, logder, particles, dims, walker_n, steps, thermal, print_interval, st):
  psi = partial(psi_alpha, alpha=alpha)
  e_l = partial(e_l_alpha, alpha=alpha)
  walkers = [cw.ConfigWalker(particles,dims) for i in range(0,walker_n)]
  samples = []
  elapsed_time = time.perf_counter() - st
  next_print = int(elapsed_time) + print_interval
  for i in range(0, steps):
    for w in walkers:
      w.metropolis_step(psi)
      if i > thermal:
        samples.append(w.current_config())
    elapsed_time = time.perf_counter() - st
    if int(elapsed_time) % print_interval == 0 and int(elapsed_time) >= next_print: 
      print(f"{elapsed_time:.1f} s: Metropolis step {i} ({i/steps * 100:.1f}%)")
      next_print = int(elapsed_time) + print_interval

  elapsed_time = time.perf_counter() - st
  print(f"{elapsed_time:.1f} s: calculating energy and stdev from {len(samples)} samples...")
  energy_last = energy
  energy = cw.sample_avg(samples, e_l)
  le_stdev = cw.sample_stdev(samples, e_l, energy)
  elapsed_time = time.perf_counter() - st
  energy_diff = np.abs(energy_last - energy)
  print(f"{elapsed_time:.1f} s: calculated energy = {energy} with stdev = {le_stdev}")
  elapsed_time = time.perf_counter() - st
  print(f"{elapsed_time:.1f} s: updating alpha...")
  alpha = update_alpha(alpha, energy, logder, e_l_alpha, samples)
  elapsed_time = time.perf_counter() - st
  print(f"{elapsed_time:.1f} s: alpha = {alpha}")
  return energy, energy_diff, le_stdev, alpha


def variational_mc(particles: int, 
                   dims: int, 
                   psi_alpha: function, 
                   e_l_alpha: function, 
                   logder: function, 
                   syst_name: str = 'unknown system'):
  """
  Run a variational Monte Carlo algorithm.\\
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
  logder: function(R, alpha)
    d/d(alpha) [ln psi_alpha], used in varying alpha\\
    Calculated analytically.
  """
  cl_args = parser.get_args()
  start_alpha = cl_args.alpha
  walker_n = cl_args.walkers
  steps = cl_args.steps
  thermal = cl_args.thermalization
  print_interval = cl_args.print
  threshold = cl_args.convergence

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
    energy, energy_diff, le_stdev, alpha = vmc_iteration(energy, alpha, e_l_alpha, psi_alpha, logder, particles, dims, walker_n, steps, thermal, print_interval, st)
    j+=1

  print(8*'=' + f" Converged with alpha = {alpha} " + 64*'=')
  print(f"Energy = {energy} with stdev = {le_stdev}")