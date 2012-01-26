import java.io.*;
import java.util.*;

/**
 * This is a driver program used to test ClusterNodes specifically
 * @author meg pirrung
 *
 */
public class TestClusterNode {
	public static void main(String[] args) {
		
		Vector sequenceData = generateSequenceData(); 
		
		Vector parseData = doParse(sequenceData);
		
		int[] mask = stripDash(parseData);
		
		System.out.println("Alignment Mask:");
		for(int i = 0; i < mask.length; i++)
			System.out.print(mask[i]);
		System.out.println();
		
		System.out.println("Alignment:");
		System.out.println(((Sequence)parseData.get(0)).first.header());
		for(int j = 0; j < parseData.size(); j++)
		{
			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					System.out.print(((String)parseData.get(j)).charAt(i));
			}		
			System.out.println();
		}
		
	}// end main

	public static Vector generateSequenceData()
	{
		Vector sData = new Vector();
		sData.add(new Sequence("Guide","  "));
		Sequence S = new Sequence("Synthetic organism","");
		S.addNodeData("  ");
		
		for (int j = 0; j < 10; j++)
		{
			String seq = S.generate(false);
//			System.out.println(seq);
			seq = seq.replace("(","");
			seq = seq.replace(")","");
			seq = seq.replace("[","");
			seq = seq.replace("]","");
			seq = seq.replace("<","");
			seq = seq.replace(">","");
			System.out.println(seq);
			Sequence SS = new Sequence("Synthetic organism",seq);
			sData.add(SS);
		}
	
		return sData;
	}

	public static Vector loadFasta()
	{
		Vector sData = new Vector();
		String temp="";
		String organism="";
		String letters="";
		sData.add(new Sequence("Guide","  "));
		BufferedReader rdr;
		try {
	    	rdr = new BufferedReader(new FileReader("5S_Rfam_Archaea_seed_Jesse_2_20_05.fasta"));
	    	
		    temp = rdr.readLine();
		    organism = temp;
		    temp = rdr.readLine();
		    
		    while(temp != null)
			{
				if(temp.charAt(0)=='>')
				{
					sData.add(new Sequence(organism.substring(1,organism.length()),letters));
					organism = temp;
					letters = "";
				}
				else
					letters+=temp;
				
				temp = rdr.readLine();
			}
		    
		    sData.add(new Sequence(organism,letters));
		    System.out.println("Read fasta file");
	    }
	    catch (IOException e) {
	 	System.out.println("Could not open fasta file");
	 	System.out.println(e);
	    }
		return sData;
	}
	
	public static Vector doParse(Vector sData)
	{
		String nodeFileName = "";
		Vector pData = new Vector();

	    String alpha = "abcdefghijklmnopqrstuvwxyz";
	    alpha += alpha;
	    alpha += alpha;
	    alpha += alpha;
	    alpha += alpha;
	    
	    pData.add(((Sequence)sData.elementAt(0)).first.header());
	    for (int i=1; i < 3; i++) 
	    {
	    	// make a *copy* of element i and work with it, not a pointer to it!
	    	Sequence S = (Sequence)sData.elementAt(i);             // focus on one sequence
	    	S.addNodeData(nodeFileName);

	    	S.parseSequence();                                   	      // parse this sequence
	    	
	    	System.out.print("Max Prob of " + S.organism + ": ");         // display parse info
	    	System.out.print(((InitialNode)S.first).optimalMaxLogProb);
	    	System.out.println();
	    	System.out.println(S.nucleotides);                            // list nucleotides

	    	if (i==-1)
	    		pData.add(((InitialNode)S.first).showParse(alpha));
	    	else
	    		pData.add(((InitialNode)S.first).showParse(S.nucleotides));

	    	
	    	System.out.println(((InitialNode)S.first).showParse(alpha));
	    	System.out.println(((InitialNode)S.first).showParse(S.nucleotides));
	    }
	    return pData;
	}
	
	public static int[] stripDash(Vector pData)
	{
		boolean found = false;
		int[] mask = new int[((String)pData.get(0)).length()];
		for(int i = 0; i<mask.length; i++)
			mask[i] = 0;
		for(int j = 0; j < mask.length; j++)
		{
			found = false;
			for(int i = 0; i < pData.size() && found != true; i++)
			{
				if(((String)pData.get(i)).charAt(j) != '-')
					found = true;
			}
			
			if(found == false)
			{
				mask[j] = 1;
			}
		}
		return mask;
	}// end align
}
