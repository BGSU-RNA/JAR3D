# JAR3D is Java-based Alignment of RNA using 3D structure information

## Motivation
Correctly predicting RNA 3D structures starting from sequence is a major, unsolved challenge in bioinformatics. An important sub-goal is the identification of the 3D structures of recurrent hairpin and internal RNA 3D motifs that appear as unpaired “loops” on secondary structure diagrams. Recurrent 3D motif functions include:

* Architectural roles introducing bends in helices (e.g., kink-turns) or changing helical twist (e.g., C-loops)
* Anchoring RNA tertiary interactions (e.g., GNRA loops and loop-receptors)
* Providing sites for proteins or small molecules to bind.
* Analysis of 3D structures (see, for example, the [RNA 3D Motif Atlas](https://rna.bgsu.edu/rna3dhub/motifs)) shows that different RNA sequences can form the same RNA 3D motif. Since not all sequence variants of a given motif are present in the 3D database, accurate methods are needed to predict which sequences are likely to form known 3D motifs. This is the purpose of JAR3D.

## Methodology
We extract all hairpin and internal loops from a [representative set of RNA 3D structures](https://rna.bgsu.edu/rna3dhub/nrlist) from the PDB and cluster them in geometrically similar families, with the results being organized into the [RNA 3D Motif Atlas](https://rna.bgsu.edu/rna3dhub/motifs).
For each recurrent motif, we construct a probabilistic model for sequence variability based on a hybrid Stochastic Context-Free Grammar/Markov Random Field (SCFG/MRF) method we developed.
To parameterize each model, we use all instances of the motif found in the Motif Atlas and knowledge of nucleotide interactions, especially isosteric basepairs and their substitution patterns.
Given the sequence of a hairpin or internal loop from a secondary structure as input, we score the sequence against each motif group to obtain both the alignment score against the SCFG/MRF model and the minimum interior edit distance to known 3D sequences, combining these into a Cutoff score which is calibrated to be 100 for perfect matches and greater than or equal to zero for reasonable matches.
If the Cutoff score is above 0, we infer that the new sequence can form the same 3D structure as in the motif group.

## Installation

The installation instructions and user manual can be found in this repository in [docs/howto.md](https://github.com/BGSU-RNA/JAR3D/blob/master/docs/howto.md).

## Data files and jar files

Data files for each release of the motif atlas and Java executable .jar files are available in the models folder in this repository.
They are also available at https://rna.bgsu.edu/data/jar3d/models/

## JAR3D web server

The JAR3D web server is ideal for relatively small investigations.
See https://rna.bgsu.edu/jar3d
