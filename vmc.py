import numpy as np

class ConfigSpace:
  """
  Configuration space for N particles in D dimensions, each of size 2L and divided into a number of subdivisions.
  """
  rng_seed = 1
  rng = np.random.default_rng(rng_seed)

  def __init__(self, particles: int, subdivisions: int, dim_size = 1., dims = 3):
    self.spatial_dims = dims
    self.dimension = particles * dims
    self.size = dim_size
    self.subdivisions = subdivisions

  def step_from(self, initial: np.ndarray, max_step: int | None = None):
    if initial.shape != (self.dimension,): raise ValueError("Initial position is wrong shape!")
    if max_step == None: max_step = max(1,self.subdivisions//1e3)
    step_axis = self.rng.integers(0, self.dimension)
    step_len = self.rng.integers(-max_step, max_step)
    step_vec = initial.copy()
    if step_vec[step_axis] - step_len <= 0: 
      step_vec[step_axis] = 0
    elif step_vec[step_axis] + step_len >= self.subdivisions-1: 
      step_vec[step_axis] = self.subdivisions-1
    else: step_vec[step_axis] += step_len
    return step_vec
  
  def __str__(self):
    particles = self.dimension//self.spatial_dims
    return f"{self.dimension}-dimensional config space (N = {particles}, D = {self.spatial_dims});"+\
          f" {self.subdivisions} positions per axis, axis length {self.size}"
  