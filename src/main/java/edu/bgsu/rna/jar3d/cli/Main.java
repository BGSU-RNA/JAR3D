package edu.bgsu.rna.jar3d.cli;

import edu.bgsu.rna.jar3d.Application;
import edu.bgsu.rna.jar3d.io.loaders.FastaLoader;
import edu.bgsu.rna.jar3d.io.loaders.QueryLoader;
import edu.bgsu.rna.jar3d.io.writers.CSVSaver;
import edu.bgsu.rna.jar3d.io.writers.ResultSaver;

public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String input = args[0];          // full path to fasta-formatted sequence file
		String models = args[1];         // path to a file telling which model files to use 
		String loopOutput = args[2];     // path to the output file giving loop-level results
		String sequenceOutput = args[3]; // path to the output file giving sequence-level results

		Application application = null;
		QueryLoader loader = null;
		ResultSaver saver = null;
		try {
			loader = new FastaLoader(input);
			saver = new CSVSaver(loopOutput, sequenceOutput);
			application = new Application(loader, saver);
			application.runAndSave("", models);
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
	}

}
