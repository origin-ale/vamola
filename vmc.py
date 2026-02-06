import time
import numpy as np
from functools import partial

import config_walk as cw
import parser

def update_alpha(alpha: float, energy: float, logder:function, e_l_alpha:function, samples: np.ndarray, gamma = 1.0):
  der = 2 * (cw.sample_avg(samples, lambda x: e_l_alpha(x,alpha)*logder(x)) - energy * cw.sample_avg(samples, logder))
  return alpha - gamma * der

def variational_mc(psi_alpha: function, e_l_alpha: function, logder: function, syst_name: str = 'unknown system'):
  """
  Run a variational Monte Carlo algorithm.
  
  Parameters
  --
  psi_alpha: function(R, alpha)
    The variational approximation for the wave function.\\
    This is a function over configuration space depending on one parameter alpha.
  e_l_alpha: function(R, alpha)
    The local energy of the system. Implements the hamiltonian and depends on the variational WF.\\
    This is a function over configuration space depending on one parameter alpha, calculated analytically.
  logder: function(R, alpha)
    d/dalpha (ln e_l_alpha), used in varying alpha\\
    This is a function over configuration space depending on one parameter alpha, calculated analytically.
  """
  cl_args = parser.get_args()
  walker_n = cl_args.walkers
  steps = cl_args.steps
  thermal = cl_args.thermalization
  print_interval = cl_args.print
  threshold = cl_args.convergence

  energy = 0
  le_stdev = np.inf
  alpha = 1.2
  st = time.perf_counter()
  j = 1
  next_print = 0
  print("Running VMC for " + f"{syst_name}.")
  print(f"Walkers: {walker_n}")
  print(f"Metropolis steps: {steps}\t\t Thermalization: {thermal}\t\t Info print interval: {print_interval} s")
  print(f"Standard deviation convergence threshold: {threshold}")
  while le_stdev > threshold:
    print(8*'-' + f" Iteration {j} " + 32*'-')
    psi = partial(psi_alpha, alpha=alpha)
    e_l = partial(e_l_alpha, alpha=alpha)
    walkers = [cw.ConfigWalker(1,1) for i in range(0,walker_n)]
    samples = []
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
    energy = cw.sample_avg(samples, e_l)
    le_stdev = cw.sample_stdev(samples, e_l, energy)
    elapsed_time = time.perf_counter() - st
    print(f"{elapsed_time:.1f} s: calculated energy = {energy} with stdev = {le_stdev}")
    old_alpha = alpha
    elapsed_time = time.perf_counter() - st
    print(f"{elapsed_time:.1f} s: updating alpha...")
    alpha = update_alpha(alpha, energy, logder, e_l_alpha, samples)
    elapsed_time = time.perf_counter() - st
    print(f"{elapsed_time:.1f} s: alpha = {alpha}")

    j+=1

  print(32*'=' + f" Converged with alpha = {alpha} " + 32*'=')
  print(f"Energy = {energy} with stdev = {le_stdev}")