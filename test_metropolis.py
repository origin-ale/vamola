import matplotlib.pyplot as plt
import numpy as np

import config_walk as cw

sigma = .1
print(f"Analytical stdev: {sigma/np.sqrt(2)}")

def test_wf(r: float):
  r0 = 1
  dist = r-r0
  return np.exp(-(dist**2)/(2* sigma**2))

walkers = [cw.ConfigWalker(1,1) for i in range(0,10)]
samples = []
for i in range(0, 100000):
  for w in walkers:
    w.metropolis_step(test_wf)
    if i >= 1000:
      samples.append(float(w.current_config()[0]))
print(f"Considering {len(samples)} samples...")
sample_hist = plt.hist(x=samples, bins=1000)
x = np.linspace(.5,1.5,1000)
y = len(samples)/275 * test_wf(x)**2
plt.plot(x,y)
avg = cw.sample_avg(samples, lambda x: x)
print(f"Average: {avg}")
sample_stdev = cw.sample_stdev(samples, lambda x: x, avg)
print(f"Stdev: {sample_stdev}")
plt.show()