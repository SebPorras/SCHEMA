SCHEMA-RICE
============

This is a software package for protein engineers. It uses protein structure and sequence information to aid researchers in designing protein recombination libraries.

SCHEMA-RASPP was developed originally developed in the laboratory of Frances H. Arnold at the California Institute of Technology. SCHEMA-RICE was developed by Weiliang Huang in the laboratory of Elizabeth Gillam at The University of Queensland. 

These tools can calculate SCHEMA energies of chimeric proteins and run the RASPP algorithm to find optimal library designs. 

SCHEMA-RICE modifies the original SCHEMA algorithm (Voigt et al., 2002) and takes into account the nature of the interaction, not just the proximity of the groups. 

$$
\large E_{1A1, 1A2} = \sum_{i \in 1A1} \sum_{j \in 1A2} C_{ij} P_{ij} M_{ij}
$$

Where:

$E$ = The disruption of interactions 

$C$ = The binary status of whether two residues fall into the interaction cuttoff distance. 

$P$ = The probability of interaction disruption. 

$M$ = The extended scoring/weighting of the disruption by physiochemical properties of amino acids. 


References:

Voigt, C. et al., "Protein building blocks preserved by recombination," Nature Structural Biology 9(7):553-558 (2002).

Meyer, M. et al., "Library analysis of SCHEMA-guided recombination," Protein Science 12:1686-1693 (2003).

Otey, C. et al., "Functional evolution and structural conservation in chimeric cytochromes P450: Calibrating a structure-guided approach," Chemistry & Biology 11:1-20 (2004)

Silberg, J. et al., "SCHEMA-guided protein recombination," Methods in Enzymology 388:35-42 (2004).

Endelman, J. et al., "Site-directed protein recombination as a shortest-path problem," Protein Engineering, Design & Selection 17(7):589-594 (2005).

Smith, M.A., Arnold, F.H., "Designing Libraries of Chimeric Proteins Using SCHEMA Recombination and RASPP", Methods in Molecular Biology 1179: 335-343 (2014)

# Installation 

1. Clone this repository to your computer. It is assumed you have Python 3.8 or higher installed. 

```
git clone https://github.com/SebPorras/SCHEMA.git
```

2. Move into the cloned directory 

```
cd SCHEMA
```

# Usage 

## SCHEMA-RICE  

There are essentially two steps to calculate optimal crossover points using SCHEMA-RICE. 

1. Generate a contact map for the proteins to be used in recombination (parental proteins). 

2. Find optimal crossover points using SCHEMA-RICE scoring and the RASPP algorithm.  

### Command line options

```
usage: rice.py [-h] -pdb PDB -msa MSA -xo XO [-pdbal PDBAL] [-chains CHAINS]
               [-min MIN] [-bin BIN] [-o output.txt] [-con contacts.txt]

Options:

    -h, --help          Show this help message and exit

    -pdb PDB            A PDB file from the Protein Data Bank

    -msa MSA            A multiple sequence alignment in ALN format (e.g. ClustalW)

    -xo XO              The number of crossovers

    -pdbal PDBAL        (Optional) In ALN format. If this argument is not provided,
                        then the PDB file's ID (e.g., 1G68) will be extracted, and
                        the sequence having that ID in the multiple sequence
                        alignment file will be used.

    -chains CHAINS      (Optional) The PDB chain identifers (e.g. -chain A B). Chains 'A' and 
                        '' are included by default.

    -min MIN            (Optional) The minimum fragment length (minus invariant positions), in
                        residues. Default min is 4.

    -bin BIN            (Optional) The width of each average mutation bin. Default bin is 1.

    -o output.txt       (Optional) Specify where you want your RASPP curve to be saved. If this 
                        option is not used, output will be printed to stdout. 

    -con contacts.txt   (Optional) You can provide an existing contact file you have previously 
                        created. If not specified, rice.py will generate a new file 
                        called contacts.txt. 
```          
### Example workflow 

#### Required files

There are two essential files you need to use the tool.

1. A **multiple sequence alignment (MSA)** of the parental proteins (proteins you wish to recombine)in ALN format without a header. An example of this is below. msa.txt is shown below. 

<pre class="fileContents">
# Multiple sequence alignment for P450s 

1A1             QVPKGLKNPPGPWGWPLIGHMLTLGKNPHLALSRMSQQYGDVLQIRIGSTPVVVLSGLDT
1A2             RVPKGLKSPPEPWGWPLLGHVLTLGKNPHLALSRMSQRYGDVLQIRIGSTPVLVLSRLDT
                :******.** ******:**:****************:**************:*** ***

1A1             IRQALVRQGDDFKGRPDLYTFTLISNGQSMSFSPDSGPVWAARRRLAQNGLKSFSIASDP
1A2             IRQALVRQGDDFKGRPDLYTSTLITDGQSLTFSTDSGPVWAARRRLAQNALNTFSIASDP
                ******************** ***::***::**.***************.*::*******
...
</pre>


2. A **PDB file** to generate the contact map. 

If one of your parent sequences in the MSA matches the 
sequence in your PDB file, you can now run the tool. For example: 

```
python rice.py -pdb 1A2.pdb -msa msa.txt -xo 6
```


However, if your PDB structure is not found in your MSA file, you will also need to provide an alignment between the PDB sequence and one of the sequences in your MSA. An example of this is shown in 1a1_2hi4_aln.txt: 

<pre class="fileContents">
2HI4            RVPKGLKSPPEPWGWPLLGHVLTLGKNPHLALSRMSQRYGDVLQIRIGSTPVLVLSRLDT
1A1             QVPKGLKNPPGPWGWPLIGHMLTLGKNPHLALSRMSQQYGDVLQIRIGSTPVVVLSGLDT
                :******.** ******:**:****************:**************:*** ***

2HI4            IRQALVRQGDDFKGRPDLYTSTLITDGQSLTFSTDSGPVWAARRRLAQNALNTFSIASDP
1A1             IRQALVRQGDDFKGRPDLYTFTLISNGQSMSFSPDSGPVWAARRRLAQNGLKSFSIASDP
                ******************** ***::***::**.***************.*::*******
</pre>

An example command would look like this: 

```
python rice.py -pdb 2HI4.pdb -msa msa.txt -pdbal 1a1_2hi4_aln.txt -xo 6
```

#### Crossover points 

The final mandatory argument is the -xo option which specifies the number of crossovers you 
wish to optimise for. 

#### Contact files

By default, rice.py will generate a new contact file called contacts.txt. However, by using the -con option, you can provide an existing contact file which you have previously created from running rice.py. This is generally recommended as it can speed up the runtime of the program. For example: 

```
python rice.py -pdb 2HI4.pdb -msa msa.txt -pdbal 1a1_2hi4_aln.txt -xo 6 -con contacts.txt
```

#### Saving your output

The -o options specifies where you would like your output to be saved. If you do not specify a file, the output will simply be printed to your stdout. For example: 

```
python rice.py -pdb 2HI4.pdb -msa msa.txt -pdbal 1a1_2hi4_aln.txt -xo 6 -o output.txt
```

The output from this command is shown below.

<pre class="fileContents">
# Minimum fragment length = 4
# Using bin width = 1
# Number of crossovers = 6
# RASPP took 48.64 secs
# RASPP found 720 results
# RASPP found 15 unique (&lt;E&gt;,&lt;m&gt;) points
# RASPP curve took 0.07 secs
# &lt;E&gt;	&lt;M&gt;	crossover points
7.0000	21.5000	304 314 337 373 426 437 
9.0000	25.0000	53 283 314 351 426 447 
8.0000	26.0000	38 85 304 314 337 426 
9.0000	26.5000	38 85 283 314 337 426 
10.5000	34.3906	243 304 314 337 426 437 
11.5000	36.2656	231 304 314 337 426 437 
12.0000	36.8281	228 304 314 337 426 437 
12.5000	37.9062	85 243 304 323 426 447 
13.5000	39.1562	85 231 304 323 426 447 
15.0000	39.7500	38 90 142 228 304 337 
14.0000	40.8281	53 110 228 304 332 426 
14.0000	41.9531	57 137 228 304 337 426 
15.0000	43.3438	81 124 188 228 304 426 
15.0000	43.5938	90 164 188 228 304 426 
16.5000	44.7344	94 188 243 279 351 426 
</pre>


#### Minimum fragment length

You can also specify a minimum fragment length using the -min option. This means that the distance bewteen crossover points will not be less than this value. By default, this value is set to 4. It is generally recommended to specify -min to prevent RASPP from choosing trival crossovers. 

```
python rice.py -pdb 1G68.pdb -msa lac-msa.txt -pdbal PSE4-1G68.txt -xo 6 -min 10 -o output.txt
```

#### Bin width 

Finally, when generating the RASPP output, you can specify the width of each average mutation bin if you wish. The default bin width is 1. For example: 

```
python rice.py -pdb 1G68.pdb -msa lac-msa.txt -pdbal PSE4-1G68.txt -xo 6 -min 10 -bin 2
```

rice.py is simply streamlines the operation of rasppcurve.py and schemacontacts.py originally provided in the SCHEMA Tools package. If you wish to use these tools independently, or understand the program better, refer to the SCHEMA Tools section below. 

## SCHEMA Tools 

The original authors have also meticulously documented the original Python tools in **schema-tools-doc.html**. 

If you are interested in exploring the use of these tools, the scripts have been updated to Python 3 but can be used exactly the same way as demonstrated in the original documentation. 

**However**, SCHEMA energy (E) is calculated using the SCHEMA-RICE algorithm which factors in the physio-chemical properties of interacting residues.  

The excerpt below details what these tools can be used for. 

>
> - Generate a contact map from a PDB file and an alignment of parent proteins
> - Calculate SCHEMA energy E and mutation m for chimeras, or an entire combinatorial library, using a contact map
> - Enumerate crossover points and compute average <E> and <m> for the resulting libraries
> - Find crossover points predicted to optimize folded, diverse proteins using RASPP
