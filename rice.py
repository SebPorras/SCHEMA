#!/usr/local/bin/python


import argparse, schema, pdb_reader, sys
from pathlib import Path


INVALID_FILE = 2
SEQUENCE_MISSING = 3
SMALL_INPUT = 4
MISSING_PARENT = 5
PDB_ALIGN_ERROR = 6
PROT_ALIGN_ERROR = 7
PARENT_ALIGN_ERROR = 8

# update if you wish to change cut off for contact
CONTACT_DISTANCE = 4.5


def main():
    args = setup_parser()

    parent_dict = read_parent_aln(args)
    pdb_contacts = generate_contacts(args, parent_dict)

    return 0


def generate_contacts(args, parent_dict: dict):
    """acts as a wrapper for original schemacontacts.py but streamlines
    it for the SCHEMA-RASPP pipeline."""

    chain_identifiers = format_chain_identifiers(args)

    residues = None
    with open(args.pdb, "r") as pdb_file:
        residues = pdb_reader.File().read(pdb_file)

    aligned_prot, aligned_pdb, pdb_key = align_pdb_residues(args, parent_dict)

    aligned_prot_is_valid(aligned_prot, pdb_key, parent_dict)
    aligned_pdb_is_valid(aligned_pdb, pdb_key, residues, chain_identifiers, args)

    # Align the residues with the parent protein
    try:
        residues = schema.alignPDBResidues(
            residues, aligned_prot, aligned_pdb, parent_dict[pdb_key], chain_identifiers
        )
    except ValueError as ve:
        print(ve)
        exit(PARENT_ALIGN_ERROR)

    # generate the contacts here
    pdb_contacts = schema.getPDBContacts(residues, CONTACT_DISTANCE)

    # display or save the contact file
    if args.cout is not None:
        with open(args.cout, "w") as contact_file:
            schema.writeContactFile(pdb_contacts, contact_file)
    else:
        schema.writeContactFile(pdb_contacts, sys.stdout)

    return pdb_contacts


def aligned_prot_is_valid(aligned_prot, pdb_key, parent_dict):
    """
    Check to make sure the parent sequence from both alignment files matches.
    """

    if aligned_prot.replace("-", "") != parent_dict[pdb_key].replace("-", ""):
        print(
            "The PDB-aligned parent and the named parent, %s, don't match!  Aborting..."
            % (pdb_key,)
        )
        exit(PROT_ALIGN_ERROR)


def aligned_pdb_is_valid(aligned_pdb, pdb_key, residues, chain_identifiers, args):
    """Check to ensure the aligned PDB sequence matches the residue sequence pulled directly from the PDB file"""

    if aligned_pdb.replace("-", "") != pdb_reader.sequence(residues, chain_identifiers):
        print(
            "The parent-aligned PDB sequence, %s, and the PDB file sequence, chain(s) %s in %s, don't match!  Aborting..."
            % (pdb_key, chain_identifiers, args.pdb)
        )
        exit(PDB_ALIGN_ERROR)


def align_pdb_residues(args, parent_dict: dict):
    """
    Because the PDB file's residue sequence may differ from those of the parents, we
    must align the PDB residues to one parent."""

    pdb_key = find_pdb_key(args)

    aligned_prot = None
    aligned_pdb = None

    if not args.pdbal:  # Just get PDB sequence from the multiple sequence alignment
        if parent_dict.get(pdb_key) is not None:
            aligned_pdb = parent_dict[pdb_key]
            aligned_prot = parent_dict[pdb_key]
        else:
            print(
                "Could not find sequence %s in the multiple sequence alignment file %s.  Aborting..."
                % (pdb_key, args.msa)
            )
            exit(SEQUENCE_MISSING)
    else:
        # Pull information from the parent/PDB alignment file.
        # Our objective is to find the sequence with the same key in both the parent MSA file and
        # the parent/PDB alignment file.
        with open(args.pdbal, "r") as par_pdb_file:
            pdb_parent_seq_list = schema.readMultipleSequenceAlignmentFile(par_pdb_file)

        pdb_parent_seqs = dict(pdb_parent_seq_list)

        # Bail out if there are fewer than 2 sequences.
        if len(pdb_parent_seqs.keys()) < 2:
            print(
                "Only found one uniquely named sequence in the PDB/parent alignment, %s.  Aborting..."
                % list(pdb_parent_seqs.keys())[0]
            )
            exit(SMALL_INPUT)

        # Find the matching key
        pdb_key = None
        for k in list(parent_dict.keys()):
            if k in pdb_parent_seqs:
                pdb_key = k

        # Bail out if no matching key is found
        if not pdb_key:
            print(
                "Could not find parents %s in PDB/parent aligned sequences %s.  Aborting..."
                % (list(parent_dict.keys()),)
            )
            exit(MISSING_PARENT)

        aligned_prot = pdb_parent_seqs[pdb_key]
        # Remove the sequence corresponding to the pdb_key, leaving only the parent sequence.
        del pdb_parent_seqs[pdb_key]

        # Take the first remaining sequence, which should be the parent sequence.
        aligned_pdb = list(pdb_parent_seqs.values())[0]

    return aligned_prot, aligned_pdb, pdb_key


def read_parent_aln(args) -> dict:
    """
    Read the alignment file to create a list of parents.
    The parents will appear in the list in the order in which they appear in the file.
    """

    with open(args.msa, "r") as msa_f:
        parents = schema.readMultipleSequenceAlignmentFile(msa_f)

    return dict(parents)


def format_chain_identifiers(args) -> list:
    """Many PDB files include multiple chains.  The chain_identifier
    list includes those chains which correspond to the protein whose
    contacts are being evaluated. Most often, chain 'A' (in the
    case of multiple chains) or chain ' ' (only one chain)
    will be the appropriate choice.
    """

    chain_identifiers = ["A", " "]
    if args.chains is not None:
        if type(args.chains) is list:
            chain_identifiers = args.chains + [" "]
        else:
            chain_identifiers = [args.chains, " "]

    return chain_identifiers


def find_pdb_key(args):
    """Will either return as None indicating that a parent
    to PDB sequence was already given, otherwise it extracts
    the PDB key from the HEADER field of the PDB file"""

    # the alignment between the reference parent
    # (indicated by reference_parent_index) and the target protein
    # sequence in the provided PDB file.  The amino acids in
    # the aligned reference parent should correspond exactly to those
    # in the msa_file above. If you don't provide a PDB alignment file,
    # the program will assume that the ID of the PDB structure
    # contained in the HEADER field corresponds to one of the
    # sequence IDs in the MSA.
    if args.pdbal is None:
        with open(args.pdb, "r") as pdb_file_name:
            return pdb_reader.File().getIDCode(pdb_file_name)

    return None


def setup_parser():
    """use the stdlib argpase to sort CLI arguments and
    return the args."""

    parser = argparse.ArgumentParser()

    # required arsgs
    parser.add_argument(
        "-pdb",
        type=validate_file,
        required=True,
        action="store",
        help="A PDB file from the Protein Data Bank",
    )
    parser.add_argument(
        "-msa",
        type=validate_file,
        required=True,
        action="store",
        help="A multiple sequence alignment in ALN format (e.g. ClustalW)",
    )
    parser.add_argument(
        "-xo", required=True, action="store", help="The number of crossovers"
    )

    # optional arguments
    # contact args
    parser.add_argument(
        "-pdbal",
        type=validate_file,
        action="store",
        help="(Optional) In ALN format. If this argument is not provided, then the PDB file's ID (e.g., 1G68) will be extracted, and the sequence having that ID in the multiple sequence alignment file will be used",
    )
    parser.add_argument(
        "-chains",
        action="store",
        help="The PDB chain identifers (e.g. -chain A B) Chains 'A' and ' ' are included by default.",
    )
    parser.add_argument(
        "-cout",
        action="store",
        help="(Optional) Will create or overwrite a contact file (eg contacts.txt)",
    )

    # RASPP args
    parser.add_argument(
        "-min",
        action="store",
        default=4,
        help="The minimum fragment length (minus invariant positions), in residues. Default min is 4",
    )
    parser.add_argument(
        "-bin", action="store", default=1, help="The width of each average mutation bin"
    )
    parser.add_argument(
        "-o", action="store", help="Will create or overwrite a file called averages.txt"
    )

    return parser.parse_args()


def validate_file(path):
    """Check that a valid file has been provided,
    otherwise exits with error code INVALID_FILE"""

    if (file := Path(path)).is_file():
        return file
    else:
        print(f"{path} is not a valid file. Aborting...")
        exit(INVALID_FILE)


if __name__ == "__main__":
    main()
