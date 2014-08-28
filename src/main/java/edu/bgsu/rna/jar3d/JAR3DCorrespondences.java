package edu.bgsu.rna.jar3d;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;

public class JAR3DCorrespondences {

	public static void main(String[] args) {
		String fastaFileName = args[0]; 
		String modelDirPath = args[1];
		String modelName = args[2];
		int rotation = Integer.parseInt(args[3]);
		String outputFileName = args[4];

		System.out.println("JAR3DCorrespondences has some issues that need to be resolved");
		
		String corrs = JAR3DMatlab.ModelCorrespondences(fastaFileName, modelDirPath, modelName, rotation);
		
		try {
			FileWriter write = new FileWriter(outputFileName,false);
			PrintWriter print_line = new PrintWriter(write);
			print_line.print(corrs);
			print_line.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
