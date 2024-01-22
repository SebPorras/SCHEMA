#!/usr/local/bin/python


import argparse, sys, time, schema, raspp, pdb_reader
from pathlib import Path

SUCCESS = 0
INVALID_FILE = 2
SEQUENCE_MISSING = 3
SMALL_INPUT = 4
MISSING_PARENT = 5
PDB_ALIGN_ERROR = 6
PROT_ALIGN_ERROR = 7
PARENT_ALIGN_ERROR = 8
FAILED_COLLAPSE_ERROR = 9

# update if you wish to change cut off for contact
CONTACT_DISTANCE = 4.5


def main():
	args = setup_parser()

	parent_list, parent_dict = read_parent_aln(args)

	parents = [p for (k,p) in parent_list]

	contacts_path = args.con
	if contacts_path is None:
		# generate new contacts
		generate_contacts(args, parent_dict)
		contacts_path = "contacts.txt"
	
	# otherwise use provided contact file
	pdb_contacts = None
	with open(contacts_path, "r") as file:
		pdb_contacts = schema.readContactFile(file)

	# sort where output will be sent 
	if args.o is None:
		output_file = sys.stdout
	else:
		output_file = open(args.o, "w")

	output_file.write(f"# Minimum fragment length = {args.min}\n")
	output_file.write(f"# Using bin width = {args.bin}\n")
	output_file.write(f"# Number of crossovers = {args.xo}\n")

	# Runs the RASPP algorithm
	find_optimal_crossovers(args, parents, pdb_contacts, output_file)
	
	if args.o is not None:
		output_file.close()

	return SUCCESS

def find_optimal_crossovers(args, parents, pdb_contacts, output_file):
	"""Wrapper for the RASPP protocols"""

	# Get the number of fragments -- one more than the number of crossovers.
	num_fragments = args.xo + 1

	# Make libraries consistent with raspp
	(new_parents, identical_sites) = raspp.collapse_parents(parents)
	new_parents_are_valid(new_parents, num_fragments, parents, args) # will exit if fails
	
	contacts = schema.getSCHEMAContacts(pdb_contacts, parents)
	energies = raspp.make_4d_energies(contacts, parents)
	
	avg_energies = raspp.calc_average_energies(energies, parents)

	tstart = time.time()
	res = raspp.RASPP(avg_energies, parents, num_fragments - 1, args.min)
	 
	output_file.write("# RASPP took %1.2f secs\n" % (time.time() - tstart,))
	output_file.write("# RASPP found %d results\n" % (len(res),))

	tstart = time.time()

	curve = raspp.curve(res, parents, args.bin)
	output_file.write("# RASPP found %d unique (<E>,<m>) points\n" % (len(curve),))
	output_file.write("# RASPP curve took %1.2f secs\n" % (time.time() - tstart,))
	output_file.write("# <E>\t<m>\tcrossover points\n")

	for average_E, average_m, crossovers in curve:
		xover_pat = "%d " * len(crossovers)
		xover_str = xover_pat % tuple(crossovers)
		output_file.write("%1.4f\t%1.4f\t%s\n" % (average_E, average_m, xover_str))


def new_parents_are_valid(new_parents, num_fragments, parents, args):
	if len(new_parents[0]) < num_fragments * args.min:
		error_msg = (
			"Minimum fragment length of %d is too large.\n%d "
			+ "fragments with length %d cannot be found in a "
			+ "sequence of length %d (with identities removed).  Aborting..."
		)
		print(error_msg % (args.min, num_fragments, args.min, len(parents[0])))
		exit(FAILED_COLLAPSE_ERROR)


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
	path = "contacts.txt"
	if args.con is not None:
		path = args.con
	
	with open(path, "w") as contact_file:
		schema.writeContactFile(pdb_contacts, contact_file)

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


def read_parent_aln(args):
	"""
	Read the alignment file to create a list of parents.
	The parents will appear in the list in the order in which they appear in the file.
	"""

	with open(args.msa, "r") as msa_f:
		parents = schema.readMultipleSequenceAlignmentFile(msa_f)

	parent_dict = dict(parents)
	return parents, parent_dict


def format_chain_identifiers(args):
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
		"-xo", type=int, required=True, action="store", help="The number of crossovers"
	)

	# optional arguments
	parser.add_argument(
		"-pdbal",
		type=validate_file,
		action="store",
		help="(Optional) In ALN format. If this argument is not provided, then the PDB file's ID (e.g., 1G68) will be extracted, and the sequence having that ID in the multiple sequence alignment file will be used",
	)
	parser.add_argument(
		"-chains",
		action="store",
		help="(Optional) The PDB chain identifers (e.g. -chain A B) Chains 'A' and ' ' are included by default.",
	)

	# RASPP args
	parser.add_argument(
		"-min",
		action="store",
		default=4,
		type=int,
		help="The minimum fragment length (minus invariant positions), in residues. Default min is 4",
	)
	parser.add_argument(
		"-bin",
		action="store",
		type=int,
		default=1,
		help = "(Optional) The width of each average mutation bin. Default bin is 1."
	)
	parser.add_argument(
		"-o", action="store", 
		metavar = "output.txt", 
		help= "(Optional) Specify where you want your RASPP curve to be saved. If this option is not used, output will be printed to stdout."
	)

	parser.add_argument(
		"-con",
		action="store",
		metavar="contacts.txt",
		type = validate_file,
		help = "(Optional) You can provide an existing contact file you have previously created. If not specified, rice.py will generate a new file called contacts.txt."
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
