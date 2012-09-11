package edu.bgsu.rna.jar3d.cli;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import edu.bgsu.rna.jar3d.Alignment;
import edu.bgsu.rna.jar3d.MotifGroup;
import edu.bgsu.rna.jar3d.Sequence;
import edu.bgsu.rna.jar3d.webJAR3D;
import edu.bgsu.rna.jar3d.query.FastaLoader;
import edu.bgsu.rna.jar3d.query.Loop;
import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.query.QueryLoader;
import edu.bgsu.rna.jar3d.query.QueryLoadingFailed;
import edu.bgsu.rna.jar3d.results.CSVSaver;
import edu.bgsu.rna.jar3d.results.LoopResult;
import edu.bgsu.rna.jar3d.results.ResultsSaver;
import edu.bgsu.rna.jar3d.results.SaveFailed;

public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String input = args[0];
		String models = args[1];
		String loopOutput = args[2];
		String sequenceOutput = args[3];

		QueryLoader loader = null;
		try {
			loader = new FastaLoader(input);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
		
		Query query = null;
		try {
			query = loader.load("");
			loader.cleanUp();
		} catch (QueryLoadingFailed e1) {
			e1.printStackTrace();
			System.exit(1);
		}
		
		List<LoopResult> results = runQuery(query, models);
		ResultsSaver saver = null;
		try {
			saver = new CSVSaver(loopOutput, sequenceOutput);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		try {
			saver.save(results);
			saver.cleanUp();
		} catch (SaveFailed e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public static List<LoopResult> runQuery(Query query, String models) {
		return MotifParse(models, query).get(0);
	}
	
	public static List<List<LoopResult>> MotifParse(String base, Query query) {
		List<List<LoopResult>> allResults = new ArrayList<List<LoopResult>>();
		for(Loop loop: query) {
			StringBuilder fasta = new StringBuilder();
			String folder = base + File.separator + loop.getType() + File.separator + "0.6";;
			
			for(String sequence: loop) {
				fasta.append(">\n"); // TODO use generated FASTA header.
				fasta.append(sequence);
				fasta.append("\n");
			}
			String fastaString = fasta.toString();
			List<LoopResult> results = MotifParse(loop.getId(), query, fastaString, folder, 
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
	public static List<LoopResult> MotifParse(long loopID, Query query, String QueryTxt, 
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
            results = Alignment.doILdbQuery((int)loopID, query, sData, modelNames, groupData, numSequences, 20);
           	
        } else {
        	results = new ArrayList<LoopResult>();
        }
	    return results;
	}

}
