# How to run JAR3D

## Installation

1. The executable version of JAR3D, named using the scheme `jar3d_DATE.jar`, can be downloaded from https://rna.bgsu.edu/data/jar3d/models/. Use the most current *DATE* if more than one is available.  Jar files from 2015 were tested with Java version 8.

2. The most current model files for hairpin loops (HL) and internal loops (IL) for each release of the [RNA 3D Motif Atlas](https://rna.bgsu.edu/rna3dhub/motifs) are found in directories at the same URL.

## Input file format

Use [FASTA format](https://en.wikipedia.org/wiki/FASTA_format), alternating a header line starting with `>` and a line with one sequence of a hairpin or internal loop. The break between strands in IL should be indicated by the `*` character. Include flanking Watson-Crick basepairs *AU*, *GC*, and *GU* in the sequence, since this is how hairpins and internal loops are extracted from 3D structures. The only allowed characters in each line are A, C, G, U, * and -.  Leave no blank lines.

**Examples:**

Single internal loop:

```
>sarcin 1
CUCAGUAU*AGAACCG
```

Multiple internal loops:

```
>sarcin 1
CUCAGUAU*AGAACCG
>sarcin 2
CUCAGUAC*GGAACCG
>sarcin 3
CUCAGUAC*GGAACUG
```

Save the file with a name like *SequenceFile.fasta* or *SequenceFile.txt*.

## Scoring a sequence file against motif groups

From the command line, run the *.jar* file like this:

```
java -jar JARFILE SEQUENCEFILE MODELLISTFILE LOOPOUTPUT SEQUENCEOUTPUT
```

where

* **JARFILE** is the name of the *.jar* file. At the time of publication, this was *jar3d_2014-12-11.jar*

* **SEQUENCEFILE** is the path to a fasta-formatted file containing the sequences to score

* **MODELLISTFILE** is the path to a text file listing the names of the model files to use. The file *all.txt* is
provided with JAR3D distributions, for example, `IL/1.13/lib/all.txt`. The model files themselves need to be in the same folder as the file listing the names of the models files to use. A number of other files need to be in the same folder as well. These are provided with the model files for RNA 3D Motif Atlas motif groups.

* **LOOPOUTPUT** is the path to a file where JAR3D will write output telling the overall score of all sequences against each model along with other diagnostics

* **SEQUENCEOUTPUT** is the path to a file where JAR3D will write output telling, for each sequence in the input file, the score and other diagnostics against each model

For example, issuing this command from the JAR3D folder containing the .jar file and release *IL/1.13* will score sequences from the 3D motif group *IL_85647.3* against all IL motif groups in release 1.13:

```
java -jar jar3d_2014-12-11.jar IL/1.13/lib/IL_85647.3.fasta IL/1.13/lib/all.txt IL_85647.3_loop.txt IL_85647.3_sequence.txt
```

## Explanation of JAR3D loop-level output

This is a **comma-delimited file** with a header line. Each line tells summary statistics covering all sequences in the input fasta file against each of the specified models. The header line is:

```
filename,motifId,%passedCutoff,meanCutoffScore,meanScore,medianScore,meanDeficit,medianDeficit,meanInteriorEditDistance,medianInteriorEditDistance,meanFullEditDistance,medianFullEditDistance,rotation
```

Here is an explanation of each field:

1. **filename** is the filename of the fasta file of sequences that were run against JAR3D

2. **motifId** is the name of the motif group that is being reported on this line￼

3. **%passedCutoff** is the percentage of sequences in the input file which fall into the acceptance region of the motif group on this line, ranging from 0 to 100

4. **meanCutoffScore** is the average CutoffScore against the present motif group, ranging from 100 downward

5. **meanScore** is the alignment score against this model, averaged across all sequences in the fasta file. The alignment score is the maximum log probability score returned by the CYK algorithm when a sequence is aligned to the SCFG/MRF model.

6. **medianScore** is the median alignment score against this model across all sequences

7. **meanDeficit** is the average alignment score deficit; for a single sequence, the deficit is the difference between the highest alignment score over all 3D instances of the motif and the score of the current sequence

8. **medianDeficit** is the median alignment score deficit

9. **meanInteriorEditDistance** refers to the edit distance between sequences in the provided fasta file and known instances of the motif from 3D structure; interior means that flanking cWW pairs are not included in the calculation of the edit distance, only the interior nucleotides. For each sequence in the fasta file, the minimum interior edit distance across all known instances from 3D is computed, then this is averaged over the sequences in the fasta file to give the meanInteriorEditDistance. So this could be called the meanMinimumInteriorEditDistance.

10. **medianInteriorEditDistance** is similar, but the median of the minimum edit distances is reported

11. **meanFullEditDistance** is similar, but includes the sequence of the flanking cWW basepairs in the calculation of the edit distance

12. **medianFullEditDistance** uses the median of the minimum edit distances

13. **rotation** is 0 for hairpin loops, 0 or 1 for internal loops, depending on whether the sequences matched the model better with the given strand order (rotation 0) or with the strands reversed (rotation 1)

## Explanation of JAR3D sequence-specific output

This is a **comma-delimited file** with a header line. The file contains one line for each sequence in the input fasta file and one line for each model, telling how the sequence scores against the model. The header line is:

```
filename,identifier,motifId,passedCutoff,meanCutoffScore,score,deficit,interiorEditDistance,fullEditDistance,rotation
```

Here is an explanation of each field on each line:

1. **filename** is the filename of the .fasta file of sequences that were run against JAR3D

2. **identifier** is the text on the line before each sequence, after the `>` character

3. **motifId** is the name of the model that is being reported on this line

4. **passedCutoff** is `true` or `false`, depending on whether the sequence is accepted or rejected by the motif group on this line

5. **meanCutoffScore** is the cutoff score of the sequence against the motif group on this line

6. **score** is the alignment score, the maximum log probability score returned by the CYK algorithm when the sequence is aligned to the SCFG/MRF model

7. **deficit** is the difference between the highest alignment score over all 3D instances of the motif and the score of the current sequence

8. **interiorEditDistance**  is the minimum edit distance between the non-flanking bases of the input sequence and each of the sequences in the current motif group, known from 3D

9. **fullEditDistance** is similar, but includes the flanking cWW basepairs in the calculation of edit distance

10. **rotation** is 0 for hairpin loops, 0 or 1 for internal loops, depending whether the sequences in the fasta file, as a whole, match the model better in the original strand order or with the strands reversed.

## Aligning sequences to a specific JAR3D motif group

Once a suitable motif group is identified, JAR3D can align the sequences to the SCFG/MRF model. The executable file *jar3dalign.jar* is used for this purpose:

`java -jar JARFILE SEQUENCEFILE MODELPATH MODELNAME ROTATION CORRESPONDENCEFILE`

where

* **JARFILE** is the name of the .jar file. At the time of publication, this was *jar3dalign_2014-12-10.jar*

* **SEQUENCEFILE** is the path to a fasta-formatted file containing the sequences to score

* **MODELPATH** is the path to the JAR3D models for the desired loop type and release, for example, *IL/1.13/lib*

* **MODELNAME** is the name of a JAR3D model, for example, **IL_85647.3**

* **ROTATION** indicates if the strands in the sequence file are reversed relative to their orientation in the￼model. This is one of the output columns of the loop and sequence-level output of the JAR3D scoring programs described above. This parameter should be 0 for non-rotated internal loops and all hairpin loops, and 1 for rotated internal loop groups.

* **CORRESPONDENCEFILE** is the path to a file where JAR3D will write output telling how each sequence position corresponds to the JAR3D model
For example, the following command will work when run from a directory containing release *IL/1.13* and the correct .jar file:

```
java -jar jar3dalign_2014-12-10.jar IL/1.13/lib/IL_85647.3.fasta IL/1.13/lib IL_85647.3 0 IL_85647.3_correspondences.txt
```

The CORRESPONDENCEFILE can be read by a human, but more readable output is produced by the Python program *fastatomodelalignment.py*. To get this program, get the JAR3D source code from [Github](https://github.com/BGSU-RNA/JAR3D) and install Python 2.7 on your computer. You do not need to compile Java programs or run Matlab programs. Run it this way:

```
python PYTHONFILE LIBDIRECTORY MOTIFID CORRESPONDENCEFILE OUTPUTFILE where
```

* **PYTHONFILE** is the path to the program fastatomodelalignment.py in the JAR3D/python directory

* **LIBDIRECTORY** is the path to the JAR3D models, such as IL/1.13/lib

* **MOTIFID** is the motif id, such as *IL_85647.3*

* **CORRESPONDENCEFILE** is the file produced by *jar3dalign* in the previous step

* **OUTPUTFILE** is an .html file which presents the alignment of the sequences in SEQUENCEFILE to the model of MOTIFID. First listed are the sequences in the motif group and how they actually correspond to the model. Second are the sequences in SEQUENCEFILE. For example, the following command turns the correspondences file produced above into an html file that displays the alignment (the final output is then the file IL_85647.3.html):

```
python jar3d/python/fastatomodelalignment.py Il/1.13/lib IL_85647.3 IL_85647.3_correspondences.txt IL_85647.3.html
```

## Getting and installing the JAR3D source code

The Matlab source code for generating SCFG/MRF probabilistic models, the Java code for parsing and aligning, and Python code for producing a presentation of the alignment are available through [Github](https://github.com/BGSU-RNA/JAR3D). The version used for the production of the article is **v1.1**.

### Running the Python code to show correspondences between sequences and a model

Install Python 2.7 and follow the instructions in Section F.6.

### Compiling the Java code to score sequences against SCFG/MRF models

This is not necessary for most users. Requires the Java Development Kit and [Maven](https://maven.apache.org/). The project is built using maven so this is done through the command line with:

```
cd jar3d
mvn package -P csv # to produce jar3d.jar, the JAR3D inference program
mvn package -P corr # to produce jar3dalign.jar, for aligning sequences to a model
```

Users interested in creating their own main class will need to edit the *pom.xml* file to use the desired *MainClass*. Packaging produces a jar file *target/jar3d-0.0.1-SNAPSHOT.jar* , which can be used to run JAR3D from the command line.

### Running the Matlab code to build SCFG/MRF models

This is not necessary for most users. In order to run the Matlab code to build SCFG/MRF models for each motif group yourself:

1. download and unzip the Matlab binary files for the desired release of the RNA 3D Motif Atlas from https://rna.bgsu.edu/data/jar3d/motifs/ and put them in a directory which we will call `JAR3DMOTIFDIRECTORY`.

2. Start Matlab and set the working directory to where you would like the models to be written; they will be written to subfolders with names such as *IL/1.13*.

3. Add to the path the directory *JAR3D/Matlab* from the Github distribution.

4. Edit the program *JAR3D/Matlab/pJAR3DMaster.m* to set the `MotifLibraryLocation` to `JAR3DMOTIFDIRECTORY`.

5. Set `JAR3Dpath` to the location of the Java code, ending with *JAR3D/target/classes*.

6. Set `Pythonpath` to the location of the Python code, ending with *JAR3D\python*.

7. Run pJAR3DMaster.m from Matlab.
