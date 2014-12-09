package edu.bgsu.rna.jar3d;

import java.io.*;
import java.util.*; 

/**
 * This is the Alignment class, it has methods that are used to facilitate aligning sequences
 * @author meg pirrung
 *
 */
public class JAR3DMatlab {

	/** 
	 * MotifCorrespondences parses a fasta file of sequences against a model
	 * and collects together correspondences between sequence positions and
	 * JAR3D model positions, returning them as a string
	 * @param fastaFileName
	 * @param modelFileName
	 * @return
	 */
	
	public static String ModelCorrespondences(String fastaFileName, String modelPath, String modelName, int rotation)
	{
		List<Sequence> sequenceData = Alignment.loadFasta(fastaFileName);
		
		MotifGroup group = new MotifGroup(modelPath, modelName);
		
		String type = group.loopType;
		
		if (rotation > 0)
			sequenceData = Alignment.reverse(sequenceData);
		
		sequenceData = Alignment.doParse(sequenceData, group.Model, 999, true, true);
		
		double[] modelScores = new double[sequenceData.size()];
		
		//Add up model scores for each sequence, find mean score, compare regular and reversed scores
	    for(int m = 1; m < sequenceData.size(); m++)
	    {
	    	double tempo = sequenceData.get(m).getMaxLogProbability(0);
	    	modelScores[m-1] =  tempo;
	    }
	    
	    int[] InteriorMinDist = Alignment.getMinEditDistance(group,sequenceData,type,rotation,true);
		int[] FullMinDist = Alignment.getMinEditDistance(group,sequenceData,type,rotation,false);
		
		boolean[] cutoffs = Alignment.getCutoffs(group,modelScores,InteriorMinDist);
		double[] cutoffscores = Alignment.getCutoffScores(group,modelScores,InteriorMinDist);

//		Alignment.displayAlignmentFASTA(sequenceData);
		
		String correspondences = "";
		
		for (int i = 1; i < sequenceData.size(); i++)
		{
			sequenceData.get(i).correspondences += "Sequence_"+i+" has_minimum_interior_edit_distance "+String.valueOf(InteriorMinDist[i-1])+"\n";
			sequenceData.get(i).correspondences += "Sequence_"+i+" has_minimum_full_edit_distance "+String.valueOf(FullMinDist[i-1])+"\n";
			sequenceData.get(i).correspondences += "Sequence_"+i+" has_cutoff_value "+String.valueOf(cutoffs[i-1])+"\n";
			sequenceData.get(i).correspondences += "Sequence_"+i+" has_cutoff_score "+String.valueOf(cutoffscores[i-1])+"\n";
		}
		
		for (int i = 1; i < sequenceData.size(); i++)
		{
            correspondences += ((Sequence)sequenceData.get(i)).correspondences;
		}
		
		correspondences = correspondences.replace("MMM",modelName);
		
		return correspondences;
	}

	public static List Align(String FastaFile, String ModelFile, int DNA, int range) 
	{
		List sequenceData = Alignment.loadFastaColumnsDNA(FastaFile,0,0,DNA); 

		sequenceData = Alignment.doParse(sequenceData, ModelFile, range);
		return sequenceData;
	}

	public static ParseData Align2(String UserDir, String FastaFile, String ModelFile, int DNA, int range) 
	{
		System.setProperty("user.dir",UserDir);
		Vector sequenceData = Alignment.loadFastaColumnsDNA(FastaFile,0,0,DNA); 
		ParseData PD = new ParseData();
		PD = Alignment.doParse2(sequenceData, ModelFile, range);
		return PD;
	}
	
	public static double[][] MotifTest(String UserDir, String loopType) 
	{
		String FASTAName = "";
		Vector<Sequence> sData;
        
		System.setProperty("user.dir",UserDir);

		Vector<String> modelNames = Sequence.getModelNames(loopType);  // read file listing model names

		System.out.println("Read " + modelNames.size() + " model names");
		
        double S[][] = new double[modelNames.size()][3*modelNames.size()];

        S[0][0] = 1;
        
		for (int m=0; m<modelNames.size();m++)                // loop through sets of sequences
	    {
	        FASTAName = ((String)modelNames.get(m)).replace(".txt",".fasta");
	        System.out.println("Aligning sequences from "+FASTAName);
	        sData = Alignment.loadFasta(FASTAName); 
	        double[] newscores = new double[2*modelNames.size()];
            if (loopType.equals("IL"))
            {
            	newscores = Alignment.getSortedILAlignment(sData, modelNames, 100);
            	for (int g=0; g < 2*modelNames.size(); g++)
            		S[m][g] = newscores[g];
            }
      	    }	    
	    return S;
	}

	
	/** 
	 * MotifTestGeneral reads a list of model names and sequence names separately
	 * It looks in the folder Sequences for sequences
	 * It looks in the folder Models for models
	 * @param UserDir
	 * @param loopType
	 * @return
	 */

	public static double[][] MotifTestGeneral(String UserDir, String loopType, String sequenceNameFile, String modelNameFile) 
	{
		String FASTAName = "";
		Vector<Sequence> sData;
        
		System.setProperty("user.dir",UserDir);
		
		Vector<String> modelNames = Sequence.readTextFile("Models" + File.separator + modelNameFile);  // read file listing model names
		System.out.println("Read " + modelNames.size() + " model names");
		Vector<String> sequenceNames = Sequence.readTextFile("Sequences" + File.separator + sequenceNameFile);  // read file listing model names
		System.out.println("Read " + sequenceNames.size() + " sequence names");
		
		double[][] S;
		if(loopType.equals("IL")){
			S = new double[sequenceNames.size()][2*modelNames.size()];
		}else{	//assume hairpin
			S = new double[sequenceNames.size()][modelNames.size()];
		}
        S[0][0] = 1;
        
		for (int m=0; m<sequenceNames.size();m++)                // loop through sets of sequences
	    {
	        FASTAName = ((String)sequenceNames.get(m));
	        System.out.println("Aligning sequences from "+FASTAName);
	        sData = Alignment.loadFasta(FASTAName); 
	        double[] newscores = new double[2*modelNames.size()];
            if (loopType.equals("IL"))
            {
            	newscores = Alignment.getSortedILAlignment(sData, modelNames, 100);
            	for (int g=0; g < 2*modelNames.size(); g++)
            		S[m][g] = newscores[g];
            }
            else if (loopType.equals("HL"))
            	newscores = Alignment.getSortedHLAlignment(sData, modelNames, 100);
            	for (int g=0; g < modelNames.size(); g++)
            		S[m][g] = newscores[g];
      	    }	    
	    return S;
	}
	
	public static Vector Display(Vector sData)
    {
		Vector<String> pData = new Vector<String>();
		
		for (int i = 0; i < sData.size(); i++)
		{
			
			pData.add(((Sequence)sData.get(i)).parseData);
			System.out.println((String)pData.get(i));
		}
		
		int[] mask = Alignment.stripDash(pData);
		String alnm = "";
		Vector alignmentVect = new Vector();

		for(int j = 0; j < pData.size(); j++)
		{
			alnm = "";
			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += pData.get(j).charAt(i);
			}		
			alnm += " "+((Sequence)sData.elementAt(j)).organism+" ";
			for(int x = 0; x < ((Sequence)sData.elementAt(j)).getMaxNodeLogProbabilitySize(); x++)
				alnm += ((Sequence)sData.elementAt(j)).getMaxNodeLogProbabilities(x).get(0);
			alignmentVect.add(alnm);
		}
		return alignmentVect;	

    }
    // Parse the sequences in the fasta file seq file against all models in the model directory
	public static double[][] MotifParseAll(String UserDir, String SeqFile, String loopType) 
	{
		Vector<Sequence> sData;
		
		System.setProperty("user.dir",UserDir);
		
		Vector<String> modelNames = Sequence.getModelNames(loopType);
		

	    sData = Alignment.loadFasta(SeqFile);
	    double S[][] = new double[sData.size()][2*modelNames.size()];
        if (loopType.equals("IL"))
        {
          	S = Alignment.getILScores(sData, modelNames, 100);
        }
	    return S;
	}
	
	public static double[] MotifParseSingle(String UserDir, String SeqFile, String ModFile) 
	{
		Vector<Sequence> sData;
		
		System.setProperty("user.dir",UserDir);
		
	    sData = Alignment.loadFasta(SeqFile);
	    
	    double S[] = new double[sData.size()];
        
        S = Alignment.getILScoresSingle(sData, ModFile, 100);
        // Please retain the following commented-out line for debugging purposes
        // Alignment.displayAlignmentFASTA(sData);

	    return S;
	}

	public static void WriteModelDists(String UserDir, String loopType,int distSize)
	{
		Vector<String> tinyModNames = new Vector<String>();
		
		Vector<String> modelNames = Sequence.getModelNames(loopType);
		for (int k=0; k< modelNames.size(); k++)
		{
			tinyModNames.add(modelNames.get(k).substring(0,5));
		}
		for(int i=0; i<modelNames.size(); i++)
		{
			String nodeFileName = modelNames.get(i);
			Sequence S = new Sequence("","");                        // blank sequence
			S.addNodeData(nodeFileName);	                         // only add node data once
		}
	}
	public static int[][] DisplayEditDists(String UserDir, String seqFile1, String seqFile2, String loopType){

		// TODO 2013-11-07 CLZ Below, the IL edit distances are called with rotation 0.  That might not be appropriate.
		
		System.setProperty("user.dir",UserDir);
		
		Vector<Sequence> seqData1 = Alignment.loadFasta(seqFile1);
		Vector<Sequence> seqData2 = Alignment.loadFasta(seqFile1);
		int[][] EditDistances;
		if(loopType.equals("IL")){
			EditDistances = SimpleAlign.calcILEditDistances(seqData1,seqData2,0,true);
		}else{
			EditDistances = SimpleAlign.calcHLEditDistances(seqData1,seqData2,true);
		}
		return EditDistances;
	}
	
	public static double[] getQuantilesFromFile(double[] Scores, String quantileFileName) {

		File distFile = new File(quantileFileName);

		Vector modDist = new Vector();
		Vector modVals = new Vector();

		try {
			FileReader inStream = new FileReader(distFile);
			BufferedReader in = new BufferedReader(inStream);
			String lineS;
			String ValueS;
			String DistS;
			double Value;
			double Dist;
			int BreakP;
			while((lineS = in.readLine()) != null){
				BreakP = lineS.indexOf(" ");
				ValueS = lineS.substring(0, BreakP);
				DistS = lineS.substring(BreakP+1);
				Dist = Double.parseDouble(DistS);
				modDist.add(Dist);
				Value = Double.parseDouble(ValueS);
				modVals.add(Value);
			}
			in.close();
		} catch(Exception e) {
			System.out.println("webJAR3D.getQuantilesFromFile: Error reading file " + quantileFileName + "\n" + e);
		}

		int n = Scores.length;
		double[] quantiles = new double[n];
		int DistLength = modVals.size();

		for(int i = 0; i < n; i++){
			int found = 0;
			int j = DistLength/2;
			int step = j;
			double current;
			double next;
			while(found==0){
				if(j==DistLength-1) {
					quantiles[i] = 1;
					found = 1;
				}
				else if(j==0){
					current = ((Double)modVals.get(j)).doubleValue();
					if(Scores[i] == current || Scores[i] > current){
						quantiles[i] = ((Double)modDist.get(j)).doubleValue();
						found = 1;
					}else{
						quantiles[i] = 0;
						found = 1;
					}
				}else{
					current = ((Double)modVals.get(j)).doubleValue();
					next = ((Double)modVals.get(j+1)).doubleValue();
					if(Scores[i] == current || (Scores[i] > current && Scores[i] < next)){
						quantiles[i] = ((Double)modDist.get(j)).doubleValue();
						found = 1;
					}else if(Scores[i] <current){
						step = Math.max(step/2,1);
						j = j - step;
					}else {
						step = Math.max(step/2,1);
						j = j + step;
					}
				}
			}
		}
		return quantiles;
	}

	public static double[] MotifTotalProbSingle(String UserDir, String SeqFile, String ModFile) 
	{
		System.setProperty("user.dir",UserDir);

		List<Sequence> sData;
		sData = Alignment.loadFasta(SeqFile);
		sData = Alignment.calculateTotalProbability(sData, ModFile, 0, false);
      
		double[] scores = new double[sData.size()-1];       // all scores computed

		for(int m = 1; m < sData.size(); m++)
		{
			scores[m-1] = sData.get(m).totalProbability;
		}

	    return scores;
	}

}