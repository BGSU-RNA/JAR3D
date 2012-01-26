import java.io.*;
import java.util.*; 

/**
 * This is the Alignment class, it has methods that are used to facilitate aligning sequences
 * @author meg pirrung
 *
 */
public class JAR3DMatlab {

	public static Vector Align(String UserDir, String FastaFile, String ModelFile, int numSequences, int DNA, int range) 
	{
		System.setProperty("user.dir",UserDir);
//		 System.out.println(System.getProperty("user.dir"));
		Vector sequenceData = Alignment.loadFastaColumnsDNA(FastaFile,0,0,DNA); 

		sequenceData = Alignment.doParse(sequenceData,numSequences,ModelFile,range);
		return sequenceData;
	}

	public static ParseData Align2(String UserDir, String FastaFile, String ModelFile, int numSequences, int DNA, int range) 
	{
		System.setProperty("user.dir",UserDir);
//		 System.out.println(System.getProperty("user.dir"));
		Vector sequenceData = Alignment.loadFastaColumnsDNA(FastaFile,0,0,DNA); 
		ParseData PD = new ParseData();
		PD = Alignment.doParse2(sequenceData,numSequences,ModelFile,range);
		return PD;
	}
	
	public static double[][] MotifTest(String UserDir, String loopType) 
	{
		int numSequences = 1000;                            // make sure this is larger than needed	
		String FASTAName = "";
		Vector sData;
		Vector scores = new Vector();
        
		System.setProperty("user.dir",UserDir);

		Vector modelNames = Sequence.getModelNames(loopType);  // read file listing model names

		System.out.println("Read " + modelNames.size() + " model names");
		
        double S[][] = new double[modelNames.size()][3*modelNames.size()];

        S[0][0] = 1;
        
		for (int m=0; m<modelNames.size();m++)                // loop through sets of sequences
	    {
	        FASTAName = ((String)modelNames.get(m)).replace(".txt",".fasta");
	        System.out.println("Aligning sequences from "+FASTAName);
	//        System.out.println("JAR3DMatlab.MotifTest: sequences\\"+FASTAName);
	        sData = Alignment.loadFasta(FASTAName); 
	        double[] newscores = new double[2*modelNames.size()];
            if (loopType.equals("IL"))
            {
            	newscores = Alignment.getSortedILAlignment(sData,modelNames,numSequences,100);
            	for (int g=0; g < 2*modelNames.size(); g++)
            		S[m][g] = newscores[g];
            }
//            else if (loopType.equals("HL"))
//            	newscores = Alignment.getSortedHLAlignment(sData,modelNames,numSequences,100);
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
		int numSequences = 10000;                     // make sure this is larger than needed	
		String FASTAName = "";
		Vector sData;
		Vector scores = new Vector();
        
		System.setProperty("user.dir",UserDir);
		
		Vector modelNames = Sequence.readTextFile("Models" + File.separator + modelNameFile);  // read file listing model names
		System.out.println("Read " + modelNames.size() + " model names");
		Vector sequenceNames = Sequence.readTextFile("Sequences" + File.separator + sequenceNameFile);  // read file listing model names
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
	//        System.out.println("JAR3DMatlab.MotifTest: sequences\\"+FASTAName);
	        sData = Alignment.loadFasta(FASTAName); 
	        double[] newscores = new double[2*modelNames.size()];
            if (loopType.equals("IL"))
            {
            	newscores = Alignment.getSortedILAlignment(sData,modelNames,numSequences,100);
            	for (int g=0; g < 2*modelNames.size(); g++)
            		S[m][g] = newscores[g];
            }
            else if (loopType.equals("HL"))
            	newscores = Alignment.getSortedHLAlignment(sData,modelNames,numSequences,100);
            	for (int g=0; g < modelNames.size(); g++)
            		S[m][g] = newscores[g];
      	    }	    
	    return S;
	}

	public static double[] MotifParse(String UserDir, String SeqFile) 
	{
		int numSequences = 10000;                            // make sure this is larger than needed	
		String FASTAName = "";
		Vector sData;
		Vector scores = new Vector();
        String loopType;
        double[] S;
        double[] newscores;
        
		System.setProperty("user.dir",UserDir);

		sData = Alignment.loadFasta(SeqFile);
        Sequence first = (Sequence)sData.elementAt(1);
        String firstLetters = first.letters;
        if (firstLetters.contains("*")){
        	loopType = "IL";
        }
        else{
        	loopType = "HL";
        }
        Vector modelNames = Sequence.getModelNames(loopType);
		
        if (loopType.equals("IL"))
        {
        	S = new double[2*modelNames.size()];
    	    newscores = new double[2*modelNames.size()];
        	newscores = Alignment.makeSortedILAlignmentHTML(sData,modelNames,numSequences,100,SeqFile,20);
           	for (int g=0; g < 2*modelNames.size(); g++)
       		S[g] = newscores[g];
        }else {  // if not IL assume HL
        	S = new double[modelNames.size()];
	    	newscores = new double[modelNames.size()];
           	newscores = Alignment.makeSortedHLAlignmentHTML(sData,modelNames,numSequences,100,SeqFile,20);
       		for (int g=0; g < modelNames.size(); g++)
       			S[g] = newscores[g];
        }
	    return S;
	}

    public static Vector Display(Vector sData)
    {
		Vector pData = new Vector();
		
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
//			for(int x = 0; x < ((Sequence)sData.elementAt(j)).maxLogProbs.size(); x++)
//				alnm += ((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(x)).get(0);
//			alnm += "   ";
			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += ((String)pData.get(j)).charAt(i);
			}		
			alnm += " "+((Sequence)sData.elementAt(j)).organism+" ";
			for(int x = 0; x < ((Sequence)sData.elementAt(j)).maxLogProbs.size(); x++)
				alnm += ((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(x)).get(0);
			alignmentVect.add(alnm);
		}
		return alignmentVect;	

    }
    // Parse the sequences in the fasta file seq file against all models in the model directory
	public static double[][] MotifParseAll(String UserDir, String SeqFile, String loopType) 
	{
		int numSequences = 10000;                            // make sure this is larger than needed	
		String FASTAName = "";
		Vector sData;
		Vector scores = new Vector();
		
		System.setProperty("user.dir",UserDir);
		
		Vector modelNames = Sequence.getModelNames(loopType);
		

	    sData = Alignment.loadFasta(SeqFile);
	    double S[][] = new double[sData.size()][2*modelNames.size()];
        if (loopType.equals("IL"))
        {
          	S = Alignment.getILScores(sData,modelNames,numSequences,100);
        }
//            else if (loopType.equals("HL"))
//            	newscores = Alignment.getSortedHLAlignment(sData,modelNames,numSequences,100);
	    return S;
	}
	public static double[] MotifParseSingle(String UserDir, String SeqFile, String ModFile) 
	{
		int numSequences = 100010;                            // make sure this is larger than needed	
		String FASTAName = "";
		Vector sData;
		Vector scores = new Vector();
		
		System.setProperty("user.dir",UserDir);
		
	    sData = Alignment.loadFasta(SeqFile);
	    
	    double S[] = new double[sData.size()];
        
        S = Alignment.getILScoresSingle(sData,ModFile,numSequences,100);

	    return S;
	}

	public static void WriteModelDists(String UserDir, String loopType,int distSize)
	{
		Vector sData;
		Vector tinyModNames = new Vector();
		
		Vector modelNames = Sequence.getModelNames(loopType);
		for (int k=0; k< modelNames.size(); k++)
		{
			tinyModNames.add(((String)modelNames.get(k)).substring(0,5));
		}
		for(int i=0; i<modelNames.size(); i++)
		{
			String nodeFileName = (String)modelNames.get(i);
			Sequence S = new Sequence("","");                        // blank sequence
			S.addNodeData(nodeFileName);	                         // only add node data once
		}
	}
	
}