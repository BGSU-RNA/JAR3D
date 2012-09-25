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
			// System.out.println(System.getProperty("user.dir"));
			DNA = (int)(Double.parseDouble(args[4]));
			Vector<Sequence> sequenceData = Alignment.loadFastaColumnsDNA(args[1],0,0,DNA); 
			numSequences = (int)(Double.parseDouble(args[3]));
			range        = (int)(Double.parseDouble(args[5]));
			sequenceData = Alignment.doParse(sequenceData,numSequences,args[2],range);
			Alignment.displayAlignmentFASTA(sequenceData,numSequences);
		}
		else
		{
		// choose sequence data and a model that goes with it
		// for index restrictions to work, the first sequence needs to be at least as long as
		// the sequence in the 3D structure from which the model was derived, and preferably
		// not much longer
		Vector<Sequence> sequenceData = Alignment.loadFasta("C:/cygwin/home/zirbel/JAR3D/sequences/IL_018_13_cWW-tSH-tHH-cSH-tWH-tHS-cWW.fasta");
		sequenceData = Alignment.doParse(sequenceData,numSequences,"C:/cygwin/home/zirbel/JAR3D/models/IL_018_13_cWW-tSH-tHH-cSH-tWH-tHS-cWW.txt",15);
		numSequences = 10;

		Alignment.displayAlignmentFASTA(sequenceData,numSequences);

		}
	}
	
}
