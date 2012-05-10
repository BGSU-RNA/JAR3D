package edu.bgsu.rna.jar3d;

import java.util.HashMap;
import java.util.List;
import java.util.Vector;

public class JAR3Database {

	//Overloaded MotifParse for new file system
		//Query should be either the full filename for the query fasta file or the full text
		//folder should be the folder with the data for the models, including loopType and version
		//modelType indicates which models to use, for example "bp".  Should be the prefix before the first "_" in model folder
		//structured is a boolean which indicates whether to use only structured models or all models
		public static List MotifParse(int loopID, String QueryTxt, String folder, String loopType, String modelType, boolean structured) 
		{
			int numSequences = 10000;                            // make sure this is larger than needed	
			String FASTAName = "";
			Vector sData;
			Vector scores = new Vector();
	        
	        double[] S;
	        List results;
	        
			System.setProperty("user.dir",folder);
	
			sData = Alignment.parseFastaText(QueryTxt,0,0);
	        
	        Vector modelNames = Sequence.getModelNames(folder, modelType, structured);
	        
	        HashMap<String,MotifGroup> groupData = webJAR3D.loadMotifGroups(folder, modelType);
	        if (loopType.equals("IL"))
	        {
	            results = Alignment.doILdbQuery(loopID, sData, modelNames, groupData, numSequences, 20);
	           	
	        }else {  // if not IL assume HL
	        	results = new Vector();
	        	S = new double[modelNames.size()];
	//	    	newscores = new double[modelNames.size()];
	//          newscores = Alignment.makeSortedHLAlignmentHTML(sData,modelNames,numSequences,100,SeqFile,20);
	//       	for (int g=0; g < modelNames.size(); g++)
	//       		S[g] = newscores[g];
	        }
		    return results;
		}

}
