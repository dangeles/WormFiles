# -*- coding: utf-8 -*-.


"""
A PCR wizard for use on the terminal command line.

It creates an object called PCR, fills in its values with the
appropriate concentrations and prints out the result to terminal.

A file is also created with the placeholder name pcr.txt and automatically
printed.

author: David Angeles-Albores
contact: dangeles at caltech dot edu
"""


class PCR(object):
    """
    An object that holds all the parts for a PCR mix.

    Attributes:
    -----------
    name - name of the enzyme to be used
    volume - volume per reaction
    n - number of reactions
    enzyme - volume of enzyme to be used per volume of reaction
    dntp - volume of dntps to be used per volume of rxn
    primers - volume of premixed primers to be used per volume of rxn
    dmso - volume of dmso to be used per volume of reaction
    dna - volume of dna to be used per volume of reaction
    buffer - volume of buffer to be used per volume of reaction
    water - volume of water to be used per volume of reaction
    """

    def __init__(self, name, n, volume):
        """
        The initialize function for this class.

        Parameters:
        name - name of the enzyme to be used, one of 'taq' or 'phusion'
        n - number of reactions to be made, type int or float
        volume - volume of the reaction, type int or float.
        """
        self.name = name.lower()
        self.volume = volume  # volume per reaction
        self.n = n*1.2  # multiply by 1.2 to prevent aliquot error
        self.enzyme = 0
        self.dntp = 1/50  # 1 ul 10mM dntp per ul rxn
        self.primers = 2.5/50  # 10uM primer per ul rxn
        self.dmso = 0
        self.dna = 0
        self.buffer = 0
        self.water = 0

    def add_enzyme(self):
        """A method to set the enzyme volume depending on self.name."""
        if self.name == 'taq':
            self.enzyme = 0.005  # ul enzyme per ul rxn
            self.buffer = 1/10  # buffer is 10x
        elif self.name == 'phusion':
            self.enzyme = .01  # ul enzyme per ul rxn
            self.buffer = 1/5  # buffer is 5x
        else:
            print('Sorry, I don\'t know how to make that mix yet')

    def add_dmso(self):
        """A method to add dmso if it's desired."""
        self.dmso = 1.5/50

    def add_dna(self, dna):
        """
        A method to specify the amount of dna to be added per vol. of rxn.

        Parameters:
        -----------
        dna - type float, the volume of dna to be added per volume of rxn.
        """
        if type(dna) is int:
            dna = float(dna)

        if (type(dna) is not float):
            raise ValueError('dna must be a number')
        self.dna = dna

    def make_mix(self):
        """
        A method to make the desired mix after all parameters are specified.

        Returns a set of 8 strings, each one of which contains
        a one line instruction of the reagent to be added.
        """
        totvol = self.volume*self.n
        sum_ingredients = (self.buffer + self.dntp + self.enzyme +
                           self.primers + self.dmso)
        self.water = self.volume - sum_ingredients*self.volume - self.dna

        s = 'Make the following mix (extra 20% added to prevent error):\n'
        print(s)

        s1 = 'Add {0:2.2f} premixed primers 10uM\n'.format(self.primers*totvol)
        print(s1)

        s2 = 'Add {0:2.2f} 10mM dNTPs\n'.format(self.dntp*totvol)
        print(s2)

        print('Add {0:2.2f} dmso'.format(self.dmso*totvol))
        s3 = 'Add {0:2.2f} dmso\n'.format(self.dmso*totvol)

        print('Add {0:2.2f} buffer'.format(self.buffer*totvol))
        s4 = 'Add {0:2.2f} buffer\n'.format(self.buffer*totvol)

        print('Add {0:2.2f} water'.format(self.water*self.n))
        s5 = 'Add {0:2.2f} water\n'.format(self.water*self.n)

        print('Add {0:2.2g} enzyme'.format(self.enzyme*totvol))
        s6 = 'Add {0:2.2g} enzyme'.format(self.enzyme*totvol)

        print('\nAdd {0:3.3} dna'.format(self.dna))
        s7 = '\nAdd {0:3.3} dna\n'.format(self.dna)

        return s, s1, s2, s3, s4, s5, s6, s7


if __name__ == '__main__':

    import argparse
    import time

    parser = argparse.ArgumentParser(description='PCR-Maker Wizard')
    parser.add_argument('name', help='A title for your PCR reaction')
    parser.add_argument('enzyme', help='Taq or Phusion', type=str)
    parser.add_argument('n', help='Number of reactions to make', type=int)
    parser.add_argument('-vol', '--vol', nargs='?',
                        help='volume per reaction, default 25', type=float)
    parser.add_argument('-d', '--dmso', nargs='?', help='Add DMSO, default no')
    parser.add_argument('-dn', '--dna', nargs='?',
                        help='Volume of DNA per rxn, default 1uL DNA/50uL rxn',
                        type=float)
    parser.add_argument('-p', '--print', help='Print', action="store_true")

    args = parser.parse_args()

    space = '-----------\n'
    name = args.name+'\n'

    if args.vol:
        rxnmsg = 'This is a master mix for {0} rxns of {1}µL\n'.format(args.n,
                                                                       args.vol
                                                                       )
    else:
        rxnmsg = 'This is a master mix for {0} rxns of {1}µL\n'.format(args.n,
                                                                       25)

    date = time.strftime("Date dd/mm/yy %d/%m/%Y") + '\n'

    wizard_info = 'PCR-Wizard Tool.\nAuthor: David Angeles-Albores.\n'
    print('')
    print(space)
    print(args.name.strip())
    print(space.strip())
    print('\033[1m' + date.strip() + '\033[0m')
    print(space.strip())
    print(rxnmsg.strip())
    print(space.strip())

    # calculations
    if args.vol:
        pcr = PCR(args.enzyme.lower(), args.n, args.vol)
    else:
        pcr = PCR(args.enzyme, args.n, 25)

    pcr.add_enzyme()

    if args.dmso:
        pcr.add_dmso()

    if args.dna:
        args.dna = float(args.dna)
        pcr.add_dna(args.dna)
    else:
        pcr.add_dna(1)

    # print results, store to file and print file
    s, s1, s2, s3, s4, s5, s6, s7 = pcr.make_mix()

    S = wizard_info+space+name+date+space+rxnmsg+space+s+s1+s2+s3+s4+s5+s6+s7

    with open('../output/pcr.txt', 'w') as f:
        f.write(S)

    import os

    if args.print:
        os.system("lpr -P SternbergLab ../output/pcr.txt")
