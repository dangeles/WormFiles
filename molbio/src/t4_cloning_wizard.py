"""
author: David Angeles-Albores.

contact: dangeles at caltech ot edu
"""
# -*- coding: utf-8 -*-.


class digest(object):
    """An object that holds a digestion mix."""

    def __init__(self, plasmid, insert='', e1='Fse1', e2='', vol=25):
        """An initialize function."""
        self.plasmid = plasmid
        self.insert = insert
        self.e1 = e1
        self.e2 = e2
        self.vol = vol

    def recipe(self, backbone=True):
        """A function to print out the digestion recipe."""
        buf = self.vol/10
        enzyme = self.vol/100
        plasmid = 4/25*self.vol
        m2 = """
After gel purification, add 1µL CIP to backbone directly, place at 37 for 30min
"""
        if backbone is True:
            message1 = 'plasmid {0}'.format(self.plasmid)
            message2 = m2
        else:
            message1 = 'insert {0}'.format(self.insert)
            message2 = 'Run a gel and gel/pcr purify.'

        if self.e2 is not '':
            water = self.vol - buf - 2*enzyme - plasmid
            s = """
Digestion for {0}:
10x Buffer = {1}µL
Enzyme {2} = {3}µL
Enzyme {4} = {3}µL
Plasmid = {5}µL
H20 = {6}µL
""".format(message1, buf, self.e1, enzyme, self.e2, plasmid, water, message2)
        else:
            water = self.vol - buf - enzyme - plasmid
            s = """
Digestion for plasmid {0}:
10x Buffer = {1}µL
Enzyme {2} = {3}µL
Plasmid = {5}µL
H20 = {6}µL

Place at 37degrees for 30 minutes.
""".format(message1, buf, self.e1, enzyme, plasmid, water, message2)

        print(s)

        return s


# class t4_ligation(object):
#     """
#     An object that holds all the parts for a T4 cloning mix.
#
#     Attributes:
#     -----------
#     name - name of the enzyme to be used
#     volume - volume per reaction
#     n - number of reactions
#     enzyme - volume of enzyme to be used per volume of reaction
#     dntp - volume of dntps to be used per volume of rxn
#     primers - volume of premixed primers to be used per volume of rxn
#     dmso - volume of dmso to be used per volume of reaction
#     dna - volume of dna to be used per volume of reaction
#     buffer - volume of buffer to be used per volume of reaction
#     water - volume of water to be used per volume of reaction
#     """
#
#     def __init__(self, name, n, volume=10):
#         """
#         The initialize function for this class.
#
#         Parameters:
#         name - name of the enzyme to be used, one of 'taq' or 'phusion'
#         n - number of reactions to be made, type int or float
#         volume - volume of the reaction, type int or float.
#         """
#         self.name = name.lower()
#         self.volume = volume  # volume per reaction
#         if n > 1:
#             self.n = n*1.2  # multiply by 1.2 to prevent aliquot error
#         else:
#             self.n = 1
#         self.enzyme = 0
#         self.vector = 0
#         self.insert = 0
#         self.buffer = 0
#         self.water = 0
#
#     def make_mix(self):
#         """
#         A method to make the desired mix after all parameters are specified.
#
#         Returns a set of 8 strings, each one of which contains
#         a one line instruction of the reagent to be added.
#         """
#         buf = self.n*self.vol/10
#         enzyme = self.n*self.vol/20
#         plasmid = self.n*self.vol/10
#         insert = 2*self.n*plasmid
#         water = self.n*self.vol - buf - enzyme - plasmid - insert
#         s = """
# T4 Ligation reaction for {0}:
# 10x Buffer = {1}µL
# Vector = {2}µL
# Insert = {3}µL
# Ligase = {4}µL
# H20 = {6}µL
#
# Let sit at room temp for 20+ minutes.
#
# Transform 3µL ligation rxn per 50µL DHalpha cells:
#                 3min ice
#                 30s 42
#                 3min ice
#                 Recovery in 500uL LB
#                 30min at 37 (skip if CARB)
#                 plate half
# """.format(self.name, buf, self.e1, enzyme, plasmid, water)
#
#         print(s)
#         return s


if __name__ == '__main__':

    import argparse
    import time
    import os

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

    wizard_info = 'T4 Cloning Tool.\nAuthor: David Angeles-Albores.\n'
    print('')
    print(space)
    print(args.name.strip())
    print(space.strip())
    print('\033[1m' + date.strip() + '\033[0m')
    print(space.strip())
    print(rxnmsg.strip())
    print(space.strip())

    # calculations

    S = wizard_info+space+name+date+space+rxnmsg+space+s+s1+s2+s3+s4+s5+s6+s7

    with open('../output/pcr.txt', 'w') as f:
        f.write(S)

    if args.print:
        os.system("lpr -P SternbergLab ../output/pcr.txt")
