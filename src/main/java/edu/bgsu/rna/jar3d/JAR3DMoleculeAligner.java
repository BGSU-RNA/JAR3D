package edu.bgsu.rna.jar3d;

import java.util.*;

/**
 * This is a driver program
 * @author meg pirrung
 *
 */
public class JAR3DMoleculeAligner {
	public static void main(String[] args) {

		int numSequences = 10;
		int DNA = 0;
		int range = 20;

		if (args.length>0)
		{
			System.setProperty("user.dir",args[0]);
			DNA = (int)(Double.parseDouble(args[4]));
			List<Sequence> sequenceData = Alignment.loadFastaColumnsDNA(args[1],0,0,DNA); 
			numSequences = (int)(Double.parseDouble(args[3]));
			range        = (int)(Double.parseDouble(args[5]));
			sequenceData = Alignment.doParse(sequenceData, args[2], range);
			sequenceData = sequenceData.subList(0, numSequences + 1);
			Alignment.displayAlignmentFASTA(sequenceData);
		}
		else
		{
			// choose sequence data and a model that goes with it
			// for index restrictions to work, the first sequence needs to be at least as long as
			// the sequence in the 3D structure from which the model was derived, and preferably
			// not much longer
			
			String FastaFile = "C:\\Users\\zirbel\\Documents\\JAR3DAlignment\\16S_sequences_from_1j5e_2AVY.fasta";
			String ModelFile = "C:\\Users\\zirbel\\Documents\\JAR3DAlignment\\16S_JAR3D1.txt";
			
			List<Sequence> sequenceData = Alignment.loadFasta(FastaFile);
			sequenceData = Alignment.doParse(sequenceData, ModelFile,15);
			sequenceData = sequenceData.subList(0, numSequences + 1);

			Alignment.displayAlignmentFASTA(sequenceData);
		}
	}
	
}
