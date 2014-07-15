#JAR3D is Java-based Alignment of RNA using 3D structures

##Motivation
Correctly predicting RNA 3D structures starting from sequence is a major, unsolved challenge in bioinformatics. An important sub-goal is the inference of the 3D structures of recurrent hairpin and internal RNA 3D motifs that appear as unpaired “loops” on secondary structure diagrams. Recurrent 3D motif functions include:

* Architectural roles introducing bends in helices (e.g. kink-turns) or changing helical twist (e.g. C-loops)
* Anchoring RNA tertiary interactions (e.g., GNRA loops and loop-receptors)
* Providing sites for proteins or small molecules to bind.
* Analysis of 3D structures (for example using WebFR3D) shows that different RNA sequences can form the same RNA 3D motif. Since not all sequence variants of a given motif are present in the 3D database, accurate methods are needed to predict which sequences are likely to form known 3D motifs. This is the purpose of JAR3D.

##Methodology
We extract all hairpin and internal loops from a non-redundant set of RNA 3D structures from the PDB/NDB and cluster them in geometrically similar families.
For each recurrent motif, we construct a probabilistic model for sequence variability based on a hybrid Stochastic Context-Free Grammar/Markov Random Field (SCFG/MRF) method we developed.
To parameterize each model, we use all instances of the motif found in the non-redundant dataset and RNA knowledge of nucleotide interactions, especially isosteric basepairs and their substitution patterns.
Given the sequence of a hairpin or internal loop from a secondary structure as input, each SCFG/MRF model calculates the probability that the sequence forms a given 3D motif. If the score is in the same range as sequences known to form the 3D structure, we infer that the new sequence can form the same 3D structure.

##Reliability
This approach correctly infers the 3D structures of nearly all structured internal loops when using as input sequences from 3D structures, a first validation step.  In many cases, a single sequence is enough to correctly infer the 3D structure. The probabilistic models for 3D motifs extracted from structurally conserved regions of 5S, 16S, or 23S rRNA were validated by scoring sequence variants from multiple sequence alignments that are different from those used to construct the models. While the NR dataset already includes a wide range of recurrent motifs, we expect new motifs to be found as new structures are solved. We have therefore created a pipeline to automatically identify new motifs in new structures as they are released by PDB and to create new SCFG/MRF models.
