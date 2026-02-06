import numpy as np
from functools import partial
import config_walk as cw
import time

def e_l_alpha(x:float, alpha: float):
  return alpha + x**2 * (0.5 - 2*alpha**2)

def psi_alpha(x:float, alpha:float):
  return np.exp(-alpha * x**2)

def logder(x:float):
  return -x**2

def update_alpha(alpha: float, energy: float, samples: np.ndarray, gamma = 1.0):
  der = 2 * (cw.sample_avg(samples, lambda x: e_l_alpha(x,alpha)*logder(x)) - energy * cw.sample_avg(samples, logder))
  return alpha - gamma * der

print_interval = 5 # Interval at which to print status, in seconds

if __name__=="__main__":
  energy = 0
  le_stdev = np.inf
  alpha = 1.2
  st = time.perf_counter()
  j = 1
  next_print = 0
  while le_stdev > 1e-4:
    psi = partial(psi_alpha, alpha=alpha)
    e_l = partial(e_l_alpha, alpha=alpha)
    walkers = [cw.ConfigWalker(1,1) for i in range(0,400)]
    samples = []
    for i in range(0, 30000):
      for w in walkers:
        w.metropolis_step(psi)
        if i > 4000:
          samples.append(w.current_config())
      elapsed_time = time.perf_counter() - st
      if int(elapsed_time) % print_interval == 0 and int(elapsed_time) >= next_print: 
        print(f"{elapsed_time:.1f} s: step {i} of iteration {j}, stdev = {le_stdev}")
        next_print = int(elapsed_time) + print_interval
    energy = cw.sample_avg(samples, e_l)
    le_stdev = cw.sample_stdev(samples, e_l, energy)
    old_alpha = alpha
    alpha = update_alpha(alpha, energy, samples)
    j+=1
print(32*'=' + f" Converged with alpha = {alpha} " + 32*'=')
print(f"Energy = {energy} with stdev = {le_stdev}")

