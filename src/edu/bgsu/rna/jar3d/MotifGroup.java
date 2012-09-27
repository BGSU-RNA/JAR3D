package edu.bgsu.rna.jar3d;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class MotifGroup implements java.io.Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	String loopType;
	String Sequences;		//Sequences seen in 3d structure - fasta format
	String Model;			//Model data in text format
	String Distribution;    //Distribution file in text format
    String[] Signature;
    int conserved;

    //Constructor that takes string folder, the full path to the folder with model 
    //and sequence folders, and the name of the group to load
    public MotifGroup(String folder, String modelType, String name) {
    	loopType = name.substring(0,1);
    	char fsep = File.separatorChar;
    	String modFolder = folder + fsep + modelType + "_models";
    	try{
    		String modelFile = modFolder+fsep+name+"_model.txt";
    		String distFile = modFolder+fsep+name+"_distribution.txt";
    		String dataFile = modFolder+fsep+name+"_data.txt";
    		String seqFile = folder+fsep+"sequences"+fsep+name+".fasta";
    		//Read in model
    		FileInputStream fstream = new FileInputStream(modelFile);
    		DataInputStream in = new DataInputStream(fstream);
    		BufferedReader br = new BufferedReader(new InputStreamReader(in));
    		String line = null;
    		StringBuilder stringBuilder = new StringBuilder();
    		String ls = System.getProperty("line.separator");
    		while( ( line = br.readLine() ) != null ) {
    	        stringBuilder.append( line );
    	        stringBuilder.append( ls );
    	    }
    		Model = stringBuilder.toString();
    		in.close();
    		//Read in sequences
    		fstream = new FileInputStream(seqFile);
    		in = new DataInputStream(fstream);
    		br = new BufferedReader(new InputStreamReader(in));
    		line = null;
    		stringBuilder = new StringBuilder();
    		ls = System.getProperty("line.separator");
    		while( ( line = br.readLine() ) != null ) {
    	        stringBuilder.append( line );
    	        stringBuilder.append( ls );
    	    }
    		Sequences = stringBuilder.toString();
    		in.close();
    		//Read in signature info
    		fstream = new FileInputStream(dataFile);
    		in = new DataInputStream(fstream);
    		br = new BufferedReader(new InputStreamReader(in));
    		String Signatures = br.readLine();
    		//TODO - fix to work with all loops, not just internal
    		Signature = new String[2];
    		int breakpoint = Signatures.indexOf(" ");
    		Signature[0] = Signatures.substring(0,breakpoint);
    		Signature[1] = Signatures.substring(breakpoint);
    		String conint = br.readLine();
    		conserved = Integer.parseInt(conint);
    		in.close();
    		//Read distribution information
    		fstream = new FileInputStream(distFile);
    		in = new DataInputStream(fstream);
    		br = new BufferedReader(new InputStreamReader(in));
    		line = null;
    		stringBuilder = new StringBuilder();
    		ls = System.getProperty("line.separator");
			while((line = br.readLine()) != null){
				stringBuilder.append( line );
    	        stringBuilder.append( ls );
			}
			Distribution = stringBuilder.toString();
			in.close();
    	}  catch (Exception e1) {
			e1.printStackTrace();
			System.out.println("Could not read input files, check path and name");
			System.out.println(e1.toString());
    	}
    }	
}