package edu.bgsu.rna.jar3d;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;

public class JAR3DCorrespondences {

	public static void main(String[] args) {
		String fastaFileName = args[0]; 
		String modelFileName = args[1];
		String modelName = args[2];
		String modelPath = "";                         // clearly wrong!
		int rotation = Integer.parseInt(args[3]);
		String outputFileName = args[4];

		System.out.println("JAR3DCorrespondences has some issues that need to be resolved");
		
		String corrs = JAR3DMatlab.ModelCorrespondences(fastaFileName, modelPath, modelFileName, rotation);
		
		try {
			FileWriter write = new FileWriter(outputFileName,false);
			PrintWriter print_line = new PrintWriter(write);
			corrs = corrs.replace("MMM", modelName);
			print_line.print(corrs);
			print_line.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
