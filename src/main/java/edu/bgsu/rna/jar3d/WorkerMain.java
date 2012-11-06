package edu.bgsu.rna.jar3d;

import java.sql.SQLException;

import edu.bgsu.rna.jar3d.io.loaders.QueryLoader;
import edu.bgsu.rna.jar3d.io.loaders.QueryLoadingFailed;
import edu.bgsu.rna.jar3d.io.writers.DBLoader;
import edu.bgsu.rna.jar3d.io.writers.DBResultSaver;
import edu.bgsu.rna.jar3d.io.writers.ResultSaver;
import edu.bgsu.rna.jar3d.io.writers.SaveFailed;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.results.LoopResult;

public class WorkerMain {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		if (args.length != 1) {
			System.out.println("Must specify a single argument, the queryId");
			System.exit(-1);
		}
		
		String queryId = args[0];
		QueryLoader loader = null;
		ResultSaver saver = null;
		Query query;
		
		try {
			loader = new DBLoader("joe", "hacker", "jdbc:mysql://localhost:3306/jar3d");
			saver = new DBResultSaver("joe", "hacker", "jdbc:mysql://localhost:3306/jar3d");
			query = loader.load(queryId);
			loader.cleanUp();
			
		} catch (SQLException e) {
			System.out.println("Could not connect to the database.");
			e.printStackTrace();
			System.exit(-1);
		} catch (QueryLoadingFailed e) {
			System.out.println("Could not load query: " + queryId);
			e.printStackTrace();
			System.exit(-1);
		}
		
//		for(Loop loop: query) {
//			LoopResult result;
//			boolean success = true;
//			try {
////				result = Alignment.analyze(loop);
//			} catch(Exception e) {
//				System.out.println(e.getMessage());
//				e.printStackTrace();
//				success = false;
//			}
////			saver.save(result, success);
//
//		}
		
		try {
			saver.cleanUp();
		} catch (SaveFailed e) {
			System.out.println("Could not clean up saving " + queryId);
			e.printStackTrace();
			System.exit(1);
		}
		System.exit(0); // Not needed but hey.
	}
}