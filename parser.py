from argparse import ArgumentParser

def vmc_args():
  walker_parser = ArgumentParser()

  walker_parser.add_argument('-a', '--alpha', action='store', type=float, default=0.8,
                             help="Starting value of the variational parameter alpha")
  walker_parser.add_argument('-w', '--walkers', action='store', type=int, default=400, 
                            help="Number of random walkers")
  walker_parser.add_argument('-s', '--steps', action='store', type=int, default=30000, 
                            help="Total Metropolis steps per walker")
  walker_parser.add_argument('-t', '--thermalization', action='store', type=int, default=4000,
                            help="Metropolis steps to discard for thermalization")
  walker_parser.add_argument('-p', '--print', action='store', type=float, default=1.,
                            help="Time interval (s) to print options during Metropolis walks")
  walker_parser.add_argument('-c', '--convergence', action='store', type=float, default=1e-6,
                            help="Convergence threshold for standard deviation of energy")
  
  return walker_parser.parse_args()

def test_args():
  test_parser = ArgumentParser()

  test_parser.add_argument('-s', '--systems', action='store', type=str, default="ho hy he",
                           help="Systems to run the testing for. 1+ of [ho: 1D harmonic oscillator, hy: hydrogen atom, he: helium atom]")
  return test_parser.parse_args()