#!/usr/local/bin/python


import sys, os, argparse
import schemacontacts as sc

from pathlib import Path

def main():

    args = setup_parser()

    print(args._get_args())

    return args


def setup_parser():
    """use the stdlib argpase to sort CLI arguments and 
    return the args."""

    parser = argparse.ArgumentParser()

    # required arsgs
    parser.add_argument("-pdb", type = validate_file,
                        required = True, action = "store", 
                        help="A PDB file from the Protein Data Bank")
    parser.add_argument("-msa", type = validate_file, 
                        required = True, action = "store", 
                        help="A multiple sequence alignment in ALN format (e.g. ClustalW)")
    parser.add_argument("-xo", required = True, action = "store", 
                        help = "The number of crossovers")

    # optional arguments
    # contact args
    parser.add_argument("-pdbal", action="store", help = "(Optional) In ALN format. If this argument is not provided, then the PDB file's ID (e.g., 1G68) will be extracted, and the sequence having that ID in the multiple sequence alignment file will be used")
    parser.add_argument("-chains", action="store", default=["A", ""], help = "The PDB chain identifers (e.g. -chain A B) Chains 'A' and ' ' are included by default.")
    parser.add_argument("-c", action="store_true", help = "(Optional) Will create or overwrite a file called contacts.txt")

    #RASPP args
    parser.add_argument("-min", action = "store", default = 4, help = "The minimum fragment length (minus invariant positions), in residues. Default min is 4")
    parser.add_argument("-bin", action = "store", default = 1, help = "The width of each average mutation bin")
    parser.add_argument("-o", action = "store_true", help = "Will create or overwrite a file called averages.txt")

    return parser.parse_args()

def validate_file(path):
    """Check that a valid file has been provided,
    other raises a FileNotFoundError"""

    if (file := Path(path)).is_file():
        return file
    else: 
        print(f"{path} is not a valid file")
        raise FileNotFoundError()


if __name__ == "__main__":
    main()

