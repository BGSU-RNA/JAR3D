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
    String[] Signature;
    Double[] Cutoffs;
//    int conserved;

    //Constructor that takes string folder, the full path to the folder with model 
    //and sequence folders, and the name of the group to load
    public MotifGroup(String folder, String modelType, String name) {
    	loopType = name.substring(0,1);
    	char fsep = File.separatorChar;
    	// String modFolder = folder + fsep + modelType + "_models";
    	// 2013-11-05 CLZ User must specify the directory in which models are to be found
    	String modFolder = folder;
    	try{
    		String modelFile = modFolder+fsep+name+"_model.txt";
    		String dataFile = modFolder+fsep+name+"_data.txt";
    		String seqFile = modFolder+fsep+name+".fasta";
    		String cutoffFile = modFolder+fsep+name+"_cutoffs.txt";
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
    		ls = "\n";                                   // needs to match how we tokenize later
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
    		int breakpoint = Signatures.indexOf(" ");
    		if(breakpoint != -1){   //Internal
    			Signature = new String[2];
    			Signature[0] = Signatures.substring(0,breakpoint);
    			Signature[1] = Signatures.substring(breakpoint);
    		}else {		//Hairpin
    			Signature = new String[1];
        		Signature[0] = Signatures;
    		}
//    		String conint = br.readLine();
//   		conserved = Integer.parseInt(conint);
    		in.close();
			//Read in cutoff values
			fstream = new FileInputStream(cutoffFile);
    		in = new DataInputStream(fstream);
    		br = new BufferedReader(new InputStreamReader(in));
    		line = br.readLine();
    		String[] Cutoffs_strings = line.split("\\t");
    		int n = Cutoffs_strings.length; 
    		Cutoffs = new Double[n];
    		for(int i = 0; i < n; i++){
    			if(Cutoffs_strings[i].equalsIgnoreCase("Inf")){
    				Cutoffs[i] = Double.MAX_VALUE;
    			}else if(Cutoffs_strings[i].equalsIgnoreCase("-Inf")){
    				Cutoffs[i] = Double.MIN_VALUE;
    			}else{
    				Cutoffs[i] = Double.valueOf(Cutoffs_strings[i]);
    			}
    		}
    	}  catch (Exception e1) {
			e1.printStackTrace();
			System.out.println("MotifGroup.MotifGroup:  Could not read model files, check path and name");
			System.out.println(e1.toString());
    	}
    }	
}