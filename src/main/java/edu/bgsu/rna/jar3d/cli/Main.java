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
		
		String input = args[0];
		String models = args[1];
		String loopOutput = args[2];
		String sequenceOutput = args[3];

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
