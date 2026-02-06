from argparse import ArgumentParser

def get_args():
  walker_parser = ArgumentParser()

  walker_parser.add_argument('-a', '--alpha', action='store', type=float, default=1.2,
                             help="Starting value of the variational parameter alpha")
  walker_parser.add_argument('-w', '--walkers', action='store', type=int, default=40, 
                            help="Number of random walkers")
  walker_parser.add_argument('-s', '--steps', action='store', type=int, default=3000, 
                            help="Total Metropolis steps per walker")
  walker_parser.add_argument('-t', '--thermalization', action='store', type=int, default=400,
                            help="Metropolis steps to discard for thermalization")
  walker_parser.add_argument('-p', '--print', action='store', type=float, default=1.,
                            help="Time interval (s) to print options during Metropolis walks")
  walker_parser.add_argument('-c', '--convergence', action='store', type=float, default=1e-4,
                            help="Convergence threshold for standard deviation of energy")
  
  return walker_parser.parse_args()
