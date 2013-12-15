package edu.bgsu.rna.jar3d;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import edu.bgsu.rna.jar3d.io.loaders.DBLoader;
import edu.bgsu.rna.jar3d.io.loaders.FastaLoader;
import edu.bgsu.rna.jar3d.io.loaders.QueryLoader;
import edu.bgsu.rna.jar3d.io.writers.CSVSaver;
import edu.bgsu.rna.jar3d.io.writers.DBResultSaver;
import edu.bgsu.rna.jar3d.io.writers.ResultSaver;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.results.LoopResult;

public class JAR3Database {

	public static void main(String[] args) {
		String models = args[0];
		String QueryID = args[1];
		String usrName = args[2];
		String pswd = args[3];
		String dbName = args[4];
		String dbConnection = "jdbc:mysql://localhost:3306/" + dbName;
		
		Application application = null;
		QueryLoader loader = null;
		ResultSaver saver = null;
		try {
			loader = new DBLoader(usrName,pswd,dbConnection);
			saver = new DBResultSaver(usrName,pswd,dbConnection);
			application = new Application(loader, saver);
			application.runAndSave(QueryID, models);
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}


	}

}
