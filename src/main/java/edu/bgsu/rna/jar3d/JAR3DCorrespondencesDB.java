package edu.bgsu.rna.jar3d;



import edu.bgsu.rna.jar3d.io.loaders.DBLoader;
import edu.bgsu.rna.jar3d.io.loaders.QueryLoader;
import edu.bgsu.rna.jar3d.io.writers.DBResultSaver;
import edu.bgsu.rna.jar3d.io.writers.ResultSaver;

public class JAR3DCorrespondencesDB {

	public static void main(String[] args) {
		String modelFolder = args[0];
		String modelName = args[1];
		String QueryID = args[2];
		String LoopID = args[3];
		String usrName = args[4];
		String pswd = args[5];
		String dbName = args[6];
		String dbConnection = "jdbc:mysql://localhost:3306/" + dbName;
		
		Application application = null;
		QueryLoader loader = null;
		ResultSaver saver = null;
		try {
			loader = new DBLoader(usrName,pswd,dbConnection);
			saver = new DBResultSaver(usrName,pswd,dbConnection);
			application = new Application(loader, saver);
			application.runAndSave(QueryID, LoopID, modelFolder, modelName);
			loader.cleanUp();
			saver.cleanUp();
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
}