package edu.bgsu.rna.jar3d.cli;

import java.io.IOException;
import java.util.List;

import edu.bgsu.rna.jar3d.query.FastaLoader;
import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.query.QueryLoader;
import edu.bgsu.rna.jar3d.query.QueryLoadingFailed;
import edu.bgsu.rna.jar3d.results.FastaSaver;
import edu.bgsu.rna.jar3d.results.LoopResult;
import edu.bgsu.rna.jar3d.results.ResultsSaver;
import edu.bgsu.rna.jar3d.results.SaveFailed;

public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length != 2) {
			System.err.println("Must give 3 arguments: input.fasta models output.fasta");
			System.exit(1);
		}
		
		String input = args[0];
		String models = args[1];
		String output = args[2];

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
		ResultsSaver saver = new FastaSaver(output);
		try {
			saver.save(results);
			saver.cleanUp();
		} catch (SaveFailed e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public static List<LoopResult> runQuery(Query query, String models) {
		return null;
	}

}
