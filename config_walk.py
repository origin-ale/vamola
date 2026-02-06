import numpy as np

class ConfigWalker:
  """
  Random walker representing the configuration space position of a system of N particles in d dimensions of size L
  """
  rng_seed = 1
  rng = np.random.default_rng(rng_seed)

  def __init__(self, particles = 1, dims = 3, scale = 1.):
    """
    Parameters
    -
    particles: int
      Number of particles in the system
    dims: int
      Number of dimensions for the system
    size: float
      Size of the system in each dimension
    """
    self.particles = particles
    self.dims = dims
    self.scale = scale
    self.coord_n = dims * particles
    self.positions = self.rng.uniform(0, scale, self.coord_n)
    self.max_step = 0.01 * scale

  def step_from(self):
    """Try stepping to a new configuration, but don't do it yet. Used in Markov chain implementation"""
    step_axis = self.rng.integers(0, self.coord_n)
    step_len = self.rng.uniform(-self.max_step, self.max_step)
    step_vec = np.zeros_like(self.positions)
    step_vec[step_axis] += step_len
    return self.positions + step_vec
  
  def move_to(self, destination: np.ndarray):
    """Move to a new configuration."""
    self.positions = destination

  def step(self):
    dest = self.step_from()
    self.move_to(dest)

  def metropolis_step(self, wf: function):
    r = self.current_config()
    rp = self.step_from()
    p = (wf(rp)**2)/(wf(r)**2)
    if p >= 1: self.move_to(rp)
    elif self.rng.random() < p: self.move_to(rp)

  def current_config(self):
    """Returns a copy of the current configuration"""
    return self.positions.copy()
  
  def __str__(self):
    return f"Walker with N = {self.particles}, d = {self.dims}, L = {self.scale} at {self.current_config()}"

def sample_avg(configs: list[np.ndarray], f: function):
  sum = 0.
  for c in configs:
    sum += f(c)
  return sum/len(configs)

def sample_stdev(configs: list[np.ndarray], f: function, avg: float):
  sum_dev = 0.
  for c in configs:
    sum_dev += (f(c) - avg)**2
  return sum_dev/len(configs)

def gauss_wf(r: np.ndarray):
  r0 = np.array([0,0,0])
  dist = np.linalg.norm(r - r0)
  return np.exp(-1000000*dist**2)

if __name__=="__main__":
  wa = ConfigWalker()
  wa.positions = np.array([0.01,0.01,0.01])
  for i in range(0,10000):
    print(wa.current_config())
    wa.metropolis_step(gauss_wf)