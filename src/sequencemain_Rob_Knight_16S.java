
import java.util.*;

/**
 * This is a driver program
 * @author meg pirrung
 *
 */
public class sequencemain_Rob_Knight_16S {
	public static void main(String[] args) {

		int numSequences = 1;

		// choose sequence data and a model that goes with it
		// for index restrictions to work, the first sequence needs to be at least as long as
		// the sequence in the 3D structure from which the model was derived, and preferably
		// not much longer
				
//		Vector sequenceData = Alignment.loadFastaColumns("rdp_download_47seqs.fasta",31,1803); 
		Vector sequenceData = Alignment.loadFasta("rdp_download_47seqs.fasta"); 
		
		numSequences = 30;

		sequenceData = Alignment.doParse(sequenceData,numSequences,"16S_model_from_2AVY.txt",30);
		
//		Alignment.displayAlignmentFASTA(sequenceData,numSequences);
		
		Alignment.displayAlignment(sequenceData,numSequences);
//		Vector aData = Alignment.getAlignment(parseData, sequenceData);
//		Alignment.printAlignment(aData, 90);
	}
}
