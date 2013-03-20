package edu.bgsu.rna.jar3d;

import java.util.*;

/**
 * 
 * @author Craig Zirbel
 *
 */
public class Tests {
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
			String modelName;
			String seqName;
			List<Sequence> sequenceData;
			double[] scores;
			double[] mlpscores;
			
			// -------------------------------------------------------------

			
			seqName = "C:/Users/zirbel/Documents/My Dropbox/BGSURNA/Motifs/lib/HL/1.0/bp-models/HL_All_Sequences_1.fasta";
			modelName = "C:/Users/zirbel/Documents/My Dropbox/BGSURNA/Motifs/lib/HL/1.0/bp-models/HL_00090.1_model.txt";

			scores = JAR3DMatlab.MotifParseSingle("C:/",seqName,modelName);
			
			
		}
	}
	
}
