import numpy as np

class ConfigWalker:
  """
  Random walker representing the configuration space position of a system of N particles in d dimensions of size L
  """
  rng_seed = 1
  rng = np.random.default_rng(rng_seed)

  def __init__(self, particles = 1, dims = 3, size = 1.):
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
    self.size = size
    self.coord_n = dims * particles
    self.positions = self.rng.uniform(0, size, self.coord_n)
    self.max_step = 0.001 * size

  def step_destination(self):
    """Try stepping to a new configuration, but don't do it yet. Used in Markov chain implementation"""
    step_axis = self.rng.integers(0, self.coord_n)
    step_len = self.rng.uniform(-self.max_step, self.max_step)
    step_vec = np.zeros_like(self.positions)
    step_vec[step_axis] += step_len
    return self.positions + step_vec
  
  def step_to(self, destination: np.ndarray):
    """Move to a new configuration."""
    self.positions = destination

  def step(self):
    dest = self.step_destination()
    self.step_to(dest)

  def current_config(self):
    """Returns a copy of the current configuration"""
    return self.positions.copy()
  
  def __str__(self):
    return f"Walker with N = {self.particles}, d = {self.dims}, L = {self.size} at {self.current_config()}"

wa = ConfigWalker()
for i in range(0, 25):
  dest = wa.step()

  

