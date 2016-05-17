"""
A wizard meant to help you keep track of brood plates in your notebook.

Author: David Angeles-Albores
contact: dangeles@caltech.edu
"""
import time


class brood_size(object):
    """
    An object that contains some relevant information for brood size assays.

    Attributes
    ----------
    name - name of the experiment
    experiments - number of treatments being handled
    replicates - number of replicates per treatment
    day - the number of passages these plates have been passaged previously
    """

    def __init__(self, name, experiments, replicates):
        """
        The initialize function.

        Parameters
        ----------
        name - name of the experiment in question
        experiments -  number of independent experiments (i.e. rnai strains)
        replicates - number of replicates per experiment
        """
        if type(name) is not str:
            raise ValueError('Name must be a string!')

        if type(experiments) not in [int, float]:
            raise ValueError('Experiments must be an int or a float')

        if type(experiments) is float:
            experiments = int(experiments)

        if type(replicates) not in [int, float]:
            raise ValueError('Replicates must be an int or a float')

        if type(replicates) is float:
            replicates = int(replicates)

        self.name = name
        self.experiments = experiments
        self.replicates = replicates
        self.date = time.strftime("%H:%M %d/%m/%Y (dd/mm/yy)")

    def summary(self):
        """A function to gather all the data into a summary."""
        info = 'Passaged plates for {0} on {1}.\n There were {2} lines'.format(
                self.name, self.date, self.experiments)
        info += ' and {0} replicates per line'.format(self.replicates)
        return info

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Brood Assay Wizard')
    parser.add_argument('name', help='A title for your experiment')
    parser.add_argument('experiments', help='Number of lines to transfer',
                        type=int)
    parser.add_argument('replicates', help='Number of replicates per line',
                        type=int)
    parser.add_argument('-p', '--print', help='Print', action="store_true")

    args = parser.parse_args()

    space = '-----------\n'
    name = args.name+'\n'
    wizard_info = 'Brood Assay Wizard Tool.\nAuthor: David Angeles-Albores.\n'

    brood_experiment = brood_size(args.name, args.experiments, args.replicates)
    # print('')
    # print(space)
    # print(args.name.strip())
    # print(space.strip())
    # print('\033[1m' + brood_experiment.date.strip() + '\033[0m')
    # print(space.strip())

    info = brood_experiment.summary()
    S = wizard_info+space+info

    print(S)

    with open('../output/brood.txt', 'w') as f:
        f.write(S)

    import os

    if args.print:
        os.system("lpr -P Sternberg_Lab ../output/brood.txt")
