
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
			Vector sequenceData = Alignment.loadFastaColumnsDNA(args[1],0,0,DNA); 
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
				
//		Vector sequenceData = Alignment.loadFasta("5S_Rfam_Archaea_seed_Jesse_2_20_05.fasta"); 
//		Vector sequenceData = Alignment.loadFasta("5S_archaeal_bacterial_combined.fasta"); 
//		Vector parseData = Alignment.doParse(sequenceData,1,"5S_model_from_1s72.txt",200);

//		Vector sequenceData = Alignment.loadFasta("Sarcin_IL_from_structures.fasta"); 
//		sequenceData = Alignment.doParse(sequenceData,100,"Sarcin.txt",30);
		
//		Vector sequenceData = Alignment.loadFasta("Internal_Loops_from_Structures.fasta"); 
//		Vector parseData = Alignment.doParse(sequenceData,100,"Internal_Loops.txt",30);

//		Vector sequenceData = Alignment.loadFasta("16S_sequences_from_2AVY_2J00.fasta"); 
//		Vector parseData = Alignment.doParse(sequenceData,2,"16S_model_from_2AVY.txt",20);

//		Vector sequenceData = Alignment.loadFasta("16S_NAR_Non_Redundant.fasta"); 
//		Vector sequenceData = Alignment.loadFasta("16S_eColi.fasta"); 
//		Vector parseData = Alignment.doParse(sequenceData,5,"16S_model_from_2AVY.txt",20);
		
//		Vector sequenceData = Alignment.loadFasta("16S_sequences_from_2J00_2AVY.fasta"); 
//		Vector parseData = Alignment.doParse(sequenceData,3,"16S_model_from_2J00.txt",20);

//		Vector sequenceData = Alignment.loadFasta("16S_sequences_from_1j5e_2AVY.fasta"); 
//		sequenceData = Alignment.doParse(sequenceData,numSequences,"16S_model_from_1j5e_mod_4.txt",10);

		Vector sequenceData = Alignment.loadFasta("C:/cygwin/home/zirbel/JAR3D/sequences/IL_018_13_cWW-tSH-tHH-cSH-tWH-tHS-cWW.fasta");
//		sequenceData = Alignment.doParse(sequenceData,numSequences,"IL_018_13_cWW-tSH-tHH-cSH-tWH-tHS-cWW.txt",15);
		sequenceData = Alignment.doParse(sequenceData,numSequences,"C:/cygwin/home/zirbel/JAR3D/models/IL_018_13_cWW-tSH-tHH-cSH-tWH-tHS-cWW.txt",15);

//		Vector sequenceData = Alignment.loadFasta("16S_sequences_from_1J5E_2QAN.fasta"); 
//		sequenceData = Alignment.doParse(sequenceData,numSequences,"16S_from_2AVY.txt",15);
//		sequenceData = Alignment.doParse(sequenceData,numSequences,"16S_model_from_1J5E_mod_4.txt",15);
				
	
		
//		Vector sequenceData = Alignment.loadFasta("23S_sequences_from_1s72_2aw4_2j01.fasta"); 
//		Vector parseData = Alignment.doParse(sequenceData,1,"23S_model_from_1s72.txt",12);
		
//		Vector sequenceData = Alignment.loadFasta("5S_Rfam_Archaea_seed_Jesse_2_20_05.fasta"); 
//		sequenceData = Alignment.doParse(sequenceData,numSequences,"5S_model_from_1s72.txt",200);

//		Vector sequenceData = Alignment.parseFastaText("> Nov 1 No 1 tWH-IL-tHS  1s72 U 1095 A 1261\nUUAAG*CGA-G-A\n> Nov 1 No 1 tWH-IL-tHS  2AW4 C  998 G 1157\nCUAAG*CGA-A-G\n> Nov 1 No 1 tWH-IL-tHS  2j01 C  998 G 1157\nCUAAG*CGA-A-G\n> Nov 1 No 1 tWH-IL-tHS  1s72 C 1456 G 1489\nCUAAG*CGAAAUG\n> Nov 1 No 1 tWH-IL-tHS  1s72 G 2773 A 2801\nGUAAG*CGA-A-A\n> Nov 1 No 1 tWH-IL-tHS  2j01 C 1351 G 1380\nCUAAG*CGA-A-G\n> Nov 1 No 1 tWH-IL-tHS  2j01 A 2738 G 2766\nAUAAC*GGA-A-G\n> Nov 1 No 1 tWH-IL-tHS  2j01 C 1577 G 1421\nCUAAG*CGA-U-G\n",0,0);
//		sequenceData = Alignment.doParse(sequenceData,numSequences,"Nov1-No1-Adjusted-tWH-IL-tHS.txt",200);

//		Vector sequenceData = Alignment.loadFasta("16S_sequences_from_2AVY_2J00.fasta"); 
//		Vector sequenceData = Alignment.loadFasta("16S_sequences_from_1j5e_2AVY.fasta"); 
//		Vector sequenceData = Alignment.loadFastaColumns("16S_NAR_Non_Redundant_1j5e_first.fasta",31,1803); 
		
		numSequences = 10;

//		Vector sequenceData = Alignment.loadFasta("LIB00003_IL_tSH-tWH-tHS.fasta");
//		sequenceData = Alignment.doParse(sequenceData,numSequences,"LIB00003_IL_tSH-tWH-tHS.txt",20);
		
		
		//sequenceData = Alignment.doParse(sequenceData,numSequences,"16S_model_from_1j5e.txt",20);
//		sequenceData = Alignment.doParse(sequenceData,numSequences,"16S_model_from_1j5e_mod_4.txt",30);
//		sequenceData = Alignment.doParse(sequenceData,numSequences,"16S_model_from_1j5e.txt",30);


		
		
		Alignment.displayAlignmentFASTA(sequenceData,numSequences);
		
//		Alignment.makeHTMLAlignment(sequenceData, numSequences);
//		Vector aData = Alignment.getAlignment(parseData, sequenceData);
//		Alignment.printAlignment(aData, 90);
		}
	}
	
}
