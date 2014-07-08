package edu.bgsu.rna.jar3d;

import java.util.*;

/**
 * This is a driver program
 * @author Craig Zirbel
 *
 */
public class JAR3DMotifAlignTest {
	public static void main(String[] args) {

		int rotation = 0;
		String fastaFileName = "C:\\Users\\zirbel\\Documents\\Motifs\\Variants_from_alignments\\IL\\IL_16330.1 IL_3U5H_057 Silva_LSU_eukaryal.fasta";
		String correspondenceFile = "C:\\Users\\zirbel\\Documents\\Motifs\\IL\\1.13\\AlignmentDiagnostics\\IL_16330.1_IL_3U5H_057_Silva_LSU_eukaryal_correspondences.txt";
		
		List<Sequence> sequenceData = Alignment.loadFasta(fastaFileName);
		
		if (rotation > 0)
			sequenceData = Alignment.reverse(sequenceData);
		
		sequenceData = Alignment.doParse(sequenceData, modelFileName, 15);

		Alignment.displayAlignmentFASTA(sequenceData);
		
		String correspondences = "";

		for (int i = 1; i < sequenceData.size(); i++)
		{
            correspondences += ((Sequence)sequenceData.get(i)).correspondences;
		}
		
		correspondences = correspondences.replace("MMM",modelName);
	}
}
