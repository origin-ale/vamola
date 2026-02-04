import numpy as np

class config_walker:
  """
  Random walker representing the configuration space position of a system of N particles in d dimensions of size L
  """
  rng_seed = 1
  rng = np.random.default_rng(rng_seed)
  next_id = 0

  def __init__(self, particles = 2, dims = 3, size = 1.):
    """
    Parameters
    -
    particles: int
      Number of particles in the system
    dims: int
      Number of dimensions for the system
    size: float or tuple
      Size of the system in each direction
    """
    self.coord_n = dims * particles
    self.positions = self.rng.random(self.coord_n) * size
    self.max_step = 0.001 * size
    self.id = self.next_id
    self.next_id += 1

  def step(self):
    step_axis = self.rng.integers(0, self.coord_n)
    step_len = self.rng.uniform(-self.max_step, self.max_step)
    self.positions[step_axis] += step_len

  def current_config(self):
    return self.positions
  
  def __str__(self):
    return f"Walker n. {self.id} at {self.current_config()}"