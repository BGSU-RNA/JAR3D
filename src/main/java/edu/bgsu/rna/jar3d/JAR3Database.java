package edu.bgsu.rna.jar3d;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import edu.bgsu.rna.jar3d.io.loaders.DBLoader;
import edu.bgsu.rna.jar3d.io.writers.DBResultSaver;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.results.LoopResult;

public class JAR3Database {

	public static void main(String[] args) {
		String base = args[0];
		String QueryID = args[1];
		String usrName = args[2];
		String pswd = args[3];
		String dbName = args[4];
		String dbConnection = "jdbc:mysql://localhost:3306/" + dbName;
		DBResultSaver rs = null;
		DBLoader db = null;
		List<List<LoopResult>> allResults = new ArrayList<List<LoopResult>>();

		try {
			db = new DBLoader(usrName,pswd,dbConnection);
		} catch(Exception e) {
			e.printStackTrace();
			System.err.println("Could not connect db for loading.");
			markFailure(usrName, pswd, dbConnection, QueryID);
		}

		Query query = null;

		try {
			query = db.load(QueryID);
		} catch(Exception e) {
			e.printStackTrace();
			System.err.println("Could not load query " + QueryID);
			markFailure(usrName, pswd, dbConnection, QueryID);
		} finally {
			db.cleanUp();
		}

		try {
			allResults = JAR3Database.MotifParse(base, query);
		} catch(Exception e){
			e.printStackTrace();
			System.err.println("Could not score query: " + QueryID);
			markFailure(usrName, pswd, dbConnection, QueryID);
		}

		try {
			rs = new DBResultSaver(usrName,pswd,dbConnection);
			for(List<LoopResult> results: allResults) {
				rs.save(results);
			}
		} catch(Exception e) {
			e.printStackTrace();
			System.err.println("Could not save: " + QueryID);
			markFailure(usrName, pswd, dbConnection, QueryID);
		} finally {
			rs.cleanUp();
		}

		System.exit(0);
	}

	private static void markFailure(String user, String password, String db, String queryId) {
		try {
			DBResultSaver rs = new DBResultSaver(user, password, db);
			rs.markFailure(queryId);
			rs.cleanUp();
		} catch(Exception e) {
			e.printStackTrace();
			System.err.println("Failed marking failure.");
		} finally {
			System.exit(-1);
		}
	}
	
	public static List<List<LoopResult>> MotifParse(String base, Query query) {
		List<List<LoopResult>> allResults = new ArrayList<List<LoopResult>>();
		for(Loop loop: query) {
			StringBuilder fasta = new StringBuilder();
			String folder = base + File.separator + loop.getTypeString() + File.separator + "1.0";

			for(String sequence: loop.getSequenceStrings()) {
				if (!sequence.isEmpty()) {
					fasta.append(">\n"); // TODO use generated FASTA header.
					fasta.append(sequence);
					fasta.append("\n");
				}
			}
			String fastaString = fasta.toString();
			List<LoopResult> results = MotifParse(loop.getId(), query, fastaString, folder, 
					loop.getTypeString(), "bp", query.onlyStructured(), loop);
			allResults.add(results);
		}
		return allResults;
	}

	//Overloaded MotifParse for new file system
	//Query should be either the full filename for the query fasta file or the full text
	//folder should be the folder with the data for the models, including loopType and version
	//modelType indicates which models to use, for example "bp".  Should be the prefix before the first "_" in model folder
	//structured is a boolean which indicates whether to use only structured models or all models
	public static List<LoopResult> MotifParse(long loopID, Query query, String QueryTxt, String folder,
			String loopType, String modelType, boolean structured, Loop loop) {
		Vector<Sequence> sData;
	    List<LoopResult> results;

		System.setProperty("user.dir", folder);

		sData = Alignment.parseFastaText(QueryTxt, 0, 0);
	    Vector<String> modelNames = Sequence.getModelNames(folder, modelType, structured);

	    HashMap<String,MotifGroup> groupData = webJAR3D.loadMotifGroups(folder, modelType);
	    if (loopType.equalsIgnoreCase("IL")) { 
	    	results = Alignment.doILdbQuery((int)loopID, query, sData, modelNames, groupData, 20, loop);

	    } else { 
	    	results = new ArrayList<LoopResult>();
	    }
		return results;
	}
}
