package edu.bgsu.rna.jar3d;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import edu.bgsu.rna.jar3d.query.Loop;
import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.query.DBLoader;
import edu.bgsu.rna.jar3d.results.DBResultSaver;
import edu.bgsu.rna.jar3d.results.LoopResult;

public class JAR3Database {
	
	public static void main(String[] args) {
		Query q;
		String base = args[0];
		String QueryID = args[1];
		String usrName = args[2];
		String pswd = args[3];
		String dbName = args[4];
		String dbConnection = "jdbc:mysql://localhost:3306/" + dbName;
		try{
			DBLoader db = new DBLoader(usrName,pswd,dbConnection);
			q = db.load(QueryID);
			List<List<LoopResult>> allResults = JAR3Database.MotifParse(base, q);
			DBResultSaver rs = new DBResultSaver(usrName,pswd,dbConnection);
			System.out.println(allResults);
			for(List<LoopResult> results: allResults){
//				System.out.println("Saving: " + results.get(0).loopId());
				rs.save(results);
			}

		}catch(Exception e){
			e.printStackTrace();
			System.out.println("Failed to load query");
		}		
	}
	
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
			System.out.println("Running: " + query.getId() + " " + loop.getId());
			List<LoopResult> results = MotifParse(loop.getId(), query, fastaString, folder, 
					loop.getType(), "bp", query.onlyStructured());
			System.out.println(" " + results);
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
	        	System.out.println("This is sane");
	            results = Alignment.doILdbQuery((int)loopID, query, sData, modelNames, groupData, numSequences, 20);
	           	
	        } else {
	        	System.out.println("Why am I seeing this");
	        	System.out.println(loopType);
	        	results = new ArrayList<LoopResult>();
	        }
		    return results;
		}
}
