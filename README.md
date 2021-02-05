# Phylogenetic Trees

In order to categorize good/bad mutations, blosum62 is used. If the score is negative, it’ll count as a non-beneficial mutation and increase the distance between the two sequences.
Instead of using actual matrices to store distances, a dictionary was used and a list of nodes was used to identify the nodes on which we can iterate.
There’s a debug mode that shows the updated distances and the nodes that were picked on each iteration, unchanged nodes distances aren’t shown as they remain the same.

## Input
The input is the name of the .txt file containing a protein family. Examples used:
- AGA2 (PF17366): A-agglutinin-binding subunit Aga2
- FB_lectin (PF07367): Fungal fruit body lectin
- v110 (PF01639): Viral family 110
- TAS2R (PF05296): Taste receptor protein (TAS2R)
- Protein C10 (PF14974)

Protein families can be obtained from https://pfam.xfam.org/family/browse and the sequences already come cropped and aligned. Chosen options:
- Seed, FASTA, Alphabetical, All upper case, Gaps as “-“ (dashes) 
So, other the examples above, more can be used in a similar fashion.

## Algorithms
- UPGMA
- Neighbor Joining

## Run Code
python3 phylo.py 

- Input: name of the .txt file containing a protein family in FASTA format), e.g.: PF17366
- Output: two .png files: one for UPGMA and one for NJ
