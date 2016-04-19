#!/usr/bin/env python

# Imports
from __future__ import print_function
import argparse
import itertools

# Constants
AA_CODONS = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTA', 'CTC', 'CTG', 'CTT'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT'],
    'X': ['TAG', 'TGA', 'TAA']
}

# Functions -------------------------------------------------------------------
def parse_arguments():
    '''Parse command line arguments.

    Parses all command line arguments and return parameters.
    '''
    parser = argparse.ArgumentParser(description="Determine if AA substitution\
                                     is possible by a single point mutation.")
    parser.add_argument("codon", type=check_codon, help="starting codon")
    parser.add_argument("amino_acid", type=check_aa, 
                        help="resulting amino acid (1 letter, stop codon = X)")
    parser.add_argument("-p", "--position", type=int, help="codon position",
                        choices=[1,2,3])
    parser.add_argument("-v", "--verbose", help="increase text output",
                        action="store_true")

    return parser.parse_args()


def check_codon(codon):
    '''Check codon input.

    Checks that codon is a 3 letter string and a valid codon.
    '''
    codon = codon.upper().replace("U","T")
    if not codon.isalpha() or len(codon) != 3 or \
      codon not in list(itertools.chain(*list(AA_CODONS.values()))):
        msg = "Codon must be a valid 3 letter string."
        raise argparse.ArgumentTypeError(msg)
    
    return codon


def check_aa(aa):
    '''Check amino acid input.

    Checks that amino acid is 1 letter string and a valid amino acid.
    '''
    if not aa.isalpha() or len(aa) != 1 or aa.upper() in "BJOUZ":
            msg = "Amino acid must be a valid 1 letter string."
            raise argparse.ArgumentTypeError(msg)
    
    return aa.upper()


def find_possible_subs(codon, aa, pos=None):
    '''Find nucleotide mutations for AA substitution, if possible.

    Takes a 3 letter codon string and single letter amino acid string
    and determines what point mutations, if any, will result in that
    amino acid substitution. Optionally, include a position integer
    to specify a specific codon position.

    If substitution is possible, returns a list of tuples length 3
    (codon position, old nucleotide, new nucleotide) for each
    potential mutation. If substitution not possible, returns None.
    '''
    old_codon = codon
    new_codons = AA_CODONS[aa]

    subs = []
    for new_codon in new_codons:
        if pos and new_codon != old_codon:
            if pos-1 == 0 and old_codon[1] == new_codon[1] and old_codon[2] == new_codon[2]:
                subs.append((pos, old_codon[pos-1], new_codon[pos-1]))
            elif pos-1 == 1 and old_codon[0] == new_codon[0] and old_codon[2] == new_codon[2]:
                subs.append((pos, old_codon[pos-1], new_codon[pos-1]))
            elif pos-1 == 2 and old_codon[0] == new_codon[0] and old_codon[1] == new_codon[1]:
                subs.append((pos, old_codon[pos-1], new_codon[pos-1]))
        elif new_codon != old_codon:
            if old_codon[1] == new_codon[1] and old_codon[2] == new_codon[2]:
                subs.append((1, old_codon[0], new_codon[0]))
            if old_codon[0] == new_codon[0] and old_codon[2] == new_codon[2]:
                subs.append((2, old_codon[1], new_codon[1]))
            if old_codon[0] == new_codon[0] and old_codon[1] == new_codon[1]:
                subs.append((3, old_codon[2], new_codon[2]))

    if not subs:
        subs = None

    return subs


def print_verbose(subs, args):
    '''Gives increased text output of results.

    Prints starting codon, amino acid substitution, positions looked
    at, and the possible mutations.
    '''
    print("--------------- Start Report ---------------")
    print("Start Codon: {0}   Substitute to: {1}".format(args.codon,
                                                         args.amino_acid))
    print("Position(s): {0}".format(args.position if args.position else "All"))
    print("Results:")
    if subs:
        if args.position:
            print("Position {0}:".format(args.position))
            for mut in subs:
                print("\t{0} -> {1}".format(mut[1],mut[2]))
        else:
            for i in [1, 2, 3]:
                found_mut = False
                print("Position {0}:".format(i))
                for mut in subs:
                    if mut[0] == i:
                        print("\t{0} -> {1}".format(mut[1],mut[2]))
                        found_mut = True
                if not found_mut:
                    print("\tNone")
    else:
        print("No Possible Point Mutations")
    print("---------------- End Report ----------------")


# Main ------------------------------------------------------------------------
def main():
    args = parse_arguments()

    subs = find_possible_subs(args.codon, args.amino_acid, args.position)
    
    if args.verbose:
        print_verbose(subs, args)
    else:
        if subs:
            for mut in subs:
                print("{0}\t{1}\t{2}".format(mut[0],mut[1],mut[2]))
        else:
            print("No Possible Point Mutations")


if __name__ == "__main__":
    main()
