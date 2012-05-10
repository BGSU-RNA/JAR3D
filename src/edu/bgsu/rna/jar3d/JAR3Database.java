package edu.bgsu.rna.jar3d;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import edu.bgsu.rna.jar3d.query.Loop;
import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.results.LoopResult;

public class JAR3Database {
	
	public static List<List<LoopResult>> MotifParse(String base, Query query) {
		List<List<LoopResult>> allResults = new ArrayList<List<LoopResult>>();
		for(Loop loop: query) {
			StringBuilder fasta = new StringBuilder();
			String folder = base + File.separator + loop.getType() + File.separator + "0.6";
			
			for(String sequence: loop) {
				fasta.append(">\n"); // TODO use generated FASTA header.
				fasta.append(sequence);
				fasta.append("\n");
			}
			String fastaString = fasta.toString();
			List<LoopResult> results = MotifParse(loop.getId(), fastaString, folder, 
					loop.getType(), "bp", query.onlyStructured());
			allResults.add(results);
		}
		return allResults;
	}

	//Overloaded MotifParse for new file system
		//Query should be either the full filename for the query fasta file or the full text
		//folder should be the folder with the data for the models, including loopType and version
		//modelType indicates which models to use, for example "bp".  Should be the prefix before the first "_" in model folder
		//structured is a boolean which indicates whether to use only structured models or all models
		public static List<LoopResult> MotifParse(long loopID, String QueryTxt, 
				String folder, String loopType, String modelType, 
				boolean structured) {
			int numSequences = 10000; // make sure this is larger than needed	
			Vector sData;	        
	        List<LoopResult> results;
	        
			System.setProperty("user.dir", folder);
	
			sData = Alignment.parseFastaText(QueryTxt, 0, 0);
	        
	        Vector modelNames = Sequence.getModelNames(folder, modelType, structured);
	        
	        HashMap<String,MotifGroup> groupData = webJAR3D.loadMotifGroups(folder, modelType);
	        if (loopType.equalsIgnoreCase("IL")) {
	            results = Alignment.doILdbQuery((int)loopID, sData, modelNames, groupData, numSequences, 20);
	           	
	        } else {
	        	results = new ArrayList<LoopResult>();
	        }
		    return results;
		}
}
