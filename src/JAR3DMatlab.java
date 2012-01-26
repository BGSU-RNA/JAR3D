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

	public static double[][] MotifTest(String UserDir, String loopType) 
	{
		int numSequences = 1000;                            // make sure this is larger than needed	
		String FASTAName = "";
		Vector sData;
		Vector scores = new Vector();
        
		System.setProperty("user.dir",UserDir);

		Vector modelNames = Sequence.getModelNames(loopType);

		System.out.println("Read " + modelNames.size() + " model names");
		
        double S[][] = new double[modelNames.size()][3*modelNames.size()];

        S[0][0] = 1;
        
		for (int m=0; m<modelNames.size();m++)                // loop through sets of sequences
	    {
	        FASTAName = ((String)modelNames.get(m)).replace(".txt",".fasta");
	        System.out.println("Aligning sequences from "+FASTAName);
	//        System.out.println("JAR3DMatlab.MotifTest: sequences\\"+FASTAName);
	        sData = Alignment.loadFasta(FASTAName); 
	        double[] newscores = new double[3*modelNames.size()];
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


}