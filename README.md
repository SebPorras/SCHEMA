SCHEMA-RICE
============

This is a software package for protein engineers. It uses protein structure and sequence information to aid researchers in designing protein recombination libraries.

SCHEMA was developed originally developed in the laboratory of Frances H. Arnold at the California Institute of Technology. 

These tools can calculate SCHEMA energies of chimeric proteins and run the RASPP algorithm to find optimal library designs. The package includes documentation and examples. 

SCHEMA-RICE modifies the original SCHEMA algorithm and takes into account the nature of the interaction, not just the proximity of the groups. 
This work was performed by Weiliang Huang in the labratory of Elizabeth Gillam at The University of Queensland. 

References:

Voigt, C. et al., "Protein building blocks preserved by recombination," Nature Structural Biology 9(7):553-558 (2002).
Meyer, M. et al., "Library analysis of SCHEMA-guided recombination," Protein Science 12:1686-1693 (2003).
Otey, C. et al., "Functional evolution and structural conservation in chimeric cytochromes P450: Calibrating a structure-guided approach," Chemistry & Biology 11:1-20 (2004)
Silberg, J. et al., "SCHEMA-guided protein recombination," Methods in Enzymology 388:35-42 (2004).
Endelman, J. et al., "Site-directed protein recombination as a shortest-path problem," Protein Engineering, Design & Selection 17(7):589-594 (2005).

# Installation 

1. Clone this repository to your computer. It is assumed you have Python 3.8 or higher installed. 

```
git clone https://github.com/SebPorras/SCHEMA.git
```


# Usage 

## SCHEMA-RICE  

There are essentially two steps to calculate optimal crossover points using SCHEMA-RICE. 

1. Generate a contact map for the proteins to be used in recombination (parental proteins). 

2. Find optimal crossover points using SCHEMA-RICE scoring and the RASPP algorithm.  


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
>
