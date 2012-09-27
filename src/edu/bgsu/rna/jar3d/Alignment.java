package edu.bgsu.rna.jar3d;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Formatter;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.results.ImmutableLoopResult;
import edu.bgsu.rna.jar3d.results.LoopResult;
import edu.bgsu.rna.jar3d.results.MutableSequenceResults;
import edu.bgsu.rna.jar3d.results.SequenceResult;

/**
 * This is the Alignment class, it has methods that are used to facilitate aligning sequences
 * @author meg pirrung
 *
 */
public class Alignment {

	public static List<Sequence> JAR3D(String UserDir, String FastaFile, String ModelFile, int numSequences, int DNA, int range) 
	{
		System.setProperty("user.dir",UserDir);
		List<Sequence> sequenceData = Alignment.loadFastaColumnsDNA(FastaFile,0,0,DNA); 
		sequenceData = Alignment.doParse(sequenceData,numSequences,ModelFile,range);
		Alignment.displayAlignmentFASTA(sequenceData,numSequences);	
		return sequenceData;
	}

	/**
	 * This method loads sequences from a FASTA file into a vector of Sequence objects called sData
	 * It reads all columns of the FASTA file
	 * @param fileName
	 * @param StartCol
	 * @param EndCol
	 * @return
	 */

	public static Vector<Sequence> loadFasta(String fileName)
	{
		return loadFastaColumns(fileName,0,0);
	}

	/**
	 * This method loads sequences from a FASTA file into a vector of Sequence objects called sData
	 * It allows you to specify the starting and ending columns to read
	 * @param fileName
	 * @param StartCol
	 * @param EndCol
	 * @return
	 */

	public static Vector<Sequence> loadFastaColumns(String fileName, int StartCol, int EndCol)
	{
		return loadFastaColumnsDNA(fileName,StartCol,EndCol,0);
	}

	/**
	 * This method loads sequences from a FASTA file into a vector of Sequence objects called sData
	 * It allows you to specify the starting and ending columns to read
	 * It converts DNA to RNA sequences if you set DNA = 1
	 * @param fileName
	 * @param StartCol
	 * @param EndCol
	 * @param DNA
	 * @return
	 */
	
	public static Vector<Sequence> loadFastaColumnsDNA(String fileName, int StartCol, int EndCol, int DNA)
		{
		Vector<Sequence> sData = new Vector<Sequence>();
		String temp="";
		String organism="";
		String letters="";
		sData.add(new Sequence(fileName,"")); // save space for header information to be added later as element 0 of this vector
		BufferedReader rdr;
		try {
			String curDir = System.getProperty("user.dir");

			try
			{
				rdr = new BufferedReader(new FileReader(fileName));
			}
			catch (FileNotFoundException e)
			{
				System.out.println("Reading sequence file from relative path");
				rdr = new BufferedReader(new FileReader(curDir + File.separator + "sequences" + File.separator + fileName));
			}
				
			temp = rdr.readLine();
			organism = temp;
			temp = rdr.readLine();

			while(temp != null)
			{
				if(temp.charAt(0)=='>')
				{
					if (DNA > 1)
					{
						String d = "";
						letters = Alignment.reverseString(letters);
                        for (int j = 0; j < 300-letters.length(); j++)
                        	d = d + '-';
                        letters = d + letters;
					}
					if (EndCol > 0)
						sData.add(new Sequence(organism.substring(1,organism.length()),letters.substring(StartCol,EndCol)));
					else
						sData.add(new Sequence(organism.substring(1,organism.length()),letters));
					organism = temp;
					letters = "";
				}
				else
				{
					if (DNA > 0)
					{
						temp = temp.replace("A", "w");
						temp = temp.replace("a", "w");
						temp = temp.replace("C", "x");
						temp = temp.replace("c", "x");
						temp = temp.replace("G", "y");
						temp = temp.replace("g", "y");
						temp = temp.replace("U", "z");
						temp = temp.replace("u", "z");
						temp = temp.replace("T", "z");
						temp = temp.replace("t", "z");
						temp = temp.replace("w", "U");
						temp = temp.replace("x", "G");
						temp = temp.replace("y", "C");
						temp = temp.replace("z", "A");
					}
					temp = temp.replace(".","-");
					letters+=temp;
					
				}

				temp = rdr.readLine();
			}

			if (DNA > 1)
			{
				String d = "";
				letters = Alignment.reverseString(letters);
                for (int j = 0; j < 300-letters.length(); j++)
                	d = d + '-';
                letters = d + letters;
			}

			
			if (EndCol > 0)
				sData.add(new Sequence(organism.substring(1,organism.length()),letters.substring(StartCol,EndCol)));
			else
				sData.add(new Sequence(organism.substring(2,organism.length()),letters));

			rdr.close();
		}
		catch (IOException e) {
			System.out.println("Could not open sequence file");
			System.out.println(e);
		}
		return sData;
	}

	/**
	 * This method sets up a vector of Sequence objects from a FASTA alignment supplied by a text field
	 * @param fastaText
	 * @param StartCol
	 * @param EndCol
	 * @return
	 */
	public static Vector<Sequence> parseFastaText(String fastaText, int StartCol, int EndCol)
	{
		Vector<Sequence> sData = new Vector<Sequence>();
		String temp="";
		String organism="";
		String letters="";

		sData.add(new Sequence("","")); // save space for header information to be added later as element 0 of this vector

		StringTokenizer st = new StringTokenizer(fastaText,"\n");

		organism = st.nextToken();
		
		while(st.hasMoreTokens())
			{
				temp = st.nextToken();
				if(temp.charAt(0)=='>')
				{
					if (EndCol > 0)
						sData.add(new Sequence(organism.substring(1,organism.length()),letters.substring(StartCol,EndCol)));
					else
						sData.add(new Sequence(organism.substring(1,organism.length()),letters));
					organism = temp;
					letters = "";
				}
				else
					letters+=temp;
			}

			sData.add(new Sequence(organism.substring(1,organism.length()),letters));

		return sData;
	}

	public static String reverseString(String S)
	{
		String T = "";
		for (int i=S.length()-1; i > 0; i--)
			T = T + S.charAt(i);
		return T;
	}

	/**
	 * This method reverses the strand orders in sequenceData and returns a new vector of sequence data
	 * @param numSequences
	 * @param sequenceData
	 * @return
	 */
	public static List<Sequence> reverse(int numSequences, List<Sequence> sequenceData)
	{
		List<Sequence> reverseSData = new Vector<Sequence>();

		String temp="";
		String organism="";
		String letters="";
		String second = "";
		String first  = "";
		
		reverseSData.add(new Sequence("","")); // save space for header information to be added later as element 0 of this vector
		
		for (int i=1; i < Math.min(sequenceData.size(),numSequences+1); i++)
		{
			temp = sequenceData.get(i).letters;
			organism = sequenceData.get(i).organism;
			
			StringTokenizer st = new StringTokenizer(temp,"*");
			second  = st.nextToken();
			first   = st.nextToken();
			letters = first + "*" + second;

			reverseSData.add(new Sequence(organism,letters));
			
		}
		
		return reverseSData;
	}

	/**
	 * This method initiates parsing.  It loops through all Sequence objects in sData, up to numSequences sequences.
	 * It adds node data just once, to a new Sequence object called S.  Then it adds the organism and letters from each Sequence object.
	 * Then it runs the parseSequence method of the Sequence class.
	 * It runs through all nodes to store max log probability information
	 * Finally, it calls the showParse method of the first node of the model, which accumulates a very, very wide version of the
	 * alignments in row-column format.  Many columns of the alignment then need to get removed later.
	 * @param sData
	 * @param numSequences
	 * @param nodeFileName
	 * @param range
	 * @return
	 */
	public static List<Sequence> doParse(List<Sequence> sData, int numSequences, String nodeFileName, int range)
	{
		Node current;
		Vector<Double> mProbs = new Vector<Double>();

		Sequence S = new Sequence("","");                        // blank sequence to use repeatedly 
		S.addNodeData(nodeFileName);	                         // only add node data once

		sData.get(0).parseData = ((InitialNode)S.first).header(); // add a header line

		Sequence firstS = sData.get(1);          // first sequence, which matches the model
		firstS.setNucleotides();
		firstS.setArrays();
		
		for (int i = 1; i < Math.min(sData.size(),numSequences+1); i++) 
		{
			S.organism = sData.get(i).organism;             // focus on one sequence
			S.letters  = sData.get(i).letters;
						
			S.setNucleotides();                                    // strip dashes from sequence
			S.setArrays();                                         // define cti, itc, convert letters to numbers

			// keep copies of cti and itc for the first sequence, which must match the model
			S.ctiFirst = new int[firstS.cti.length];            
			S.itcFirst = new int[firstS.itc.length];
			
			for (int j = 0; j < S.ctiFirst.length; j++)
				S.ctiFirst[j] = firstS.cti[j];
			for (int j = 0; j < S.itcFirst.length; j++)
				S.itcFirst[j] = firstS.itc[j];

			S.parseSequence(range);                             	      // parse this sequence

			mProbs = new Vector<Double>();
			current = S.first;
			while(current != null)
			{
				mProbs.add(new Double(current.optimalMaxLogProb));
				current =  current.next;
			}
			sData.get(i).appendProbabilities(mProbs);
			if (S.first.optimalMaxLogProb < -99999999)
			{
				current = S.first.next;
				while (current != null)
				{
					if (current.previous.currentMaxLogProb < -99999999 && current.currentMaxLogProb > -999999999)
					{
					}
					current = current.next;				
				}
			}
			
			sData.get(i).parseData = ((InitialNode)S.first).showParse(S.nucleotides);
			
			String correspondences = ((InitialNode)S.first).showCorrespondences(S.nucleotides);

			correspondences = correspondences.replace("JAR3D_aligns_to", "aligns_to_JAR3D");
			
			String SF = "Sequence_"+i+"_"+S.organism;
			SF = SF.replace(" ","_");
			correspondences = correspondences.replace("SSS",SF);
			
			String NF = "MMM";                           // model name goes here
			correspondences = correspondences.replace("MMM", NF);
			sData.get(i).correspondences = correspondences;
		}
		return sData;
	}

	/**
	 * This method generates synthetic sequences based on the model selected
	 * @param numSequences is the number of sequences to synthesize
	 * @param modelNum is the number corresponding to which model to use (will probably be a text file soon)
	 * @return
	 */
	public static Vector<Sequence> generateSyntheticData(int numSequences,int range)
	{
		Sequence S = new Sequence("Synthetic organism","AAAA");
		S.addNodeData("");

		Vector<Sequence> sData = new Vector<Sequence>();

		sData.add(new Sequence("",""));

		for (int i = 0; i < numSequences; i++)
		{
			String seq = S.generate(false);
			System.out.println(seq);
			seq = seq.replace("(","");
			seq = seq.replace(")","");
			seq = seq.replace("[","");
			seq = seq.replace("]","");
			seq = seq.replace("<","");
			seq = seq.replace(">","");
			seq = seq.replace("{","");
			seq = seq.replace("}","");
			seq = seq.replace("-","");
			sData.add(new Sequence("Synthetic organism",seq));
		}

		return sData;
	}


	/**
	 * This method strips out repetitive dashes used to align sequences.  This is how the very, very wide alignment gets narrowed.
	 * @param pData this is a vector of plain string sequences
	 * @return
	 */
	public static int[] stripDash(List<String> pData)
	{
		boolean found = false;
		int[] mask = new int[pData.get(0).length()];
		for(int i = 0; i<mask.length; i++)
			mask[i] = 0;
		for(int j = 0; j < mask.length; j++)                       // look at each column
		{
			found = false;                                         // has a non-dash been found?
			
			for(int i = 0; i < pData.size() && found != true; i++)
			{
				if(pData.get(i).charAt(j) != '-')
					found = true;                                  // a non-dash has been found
			}
			if(found == false)
			{
				mask[j] = 1;
			}
		}
		for(int j = 0; j < mask.length; j++)
		{
			found = false;
			for(int i = 0; i < pData.size() && found != true; i++)
			{
				if(pData.get(i).charAt(j) != '(')
					found = true;
			}

			if(found == false)
			{
				mask[j] = 1;
			}
		}
		for(int j = 0; j < mask.length; j++)
		{
			found = false;
			for(int i = 0; i < pData.size() && found != true; i++)
			{
				if(pData.get(i).charAt(j) != ')')
					found = true;
			}

			if(found == false)
			{
				mask[j] = 1;
			}
		}
		for(int j = 0; j < mask.length; j++)
		{
			found = false;
			for(int i = 0; i < pData.size() && found != true; i++)
			{
				if(pData.get(i).charAt(j) != '[')
					found = true;
			}

			if(found == false)
			{
				mask[j] = 1;
			}
		}
		for(int j = 0; j < mask.length; j++)
		{
			found = false;
			for(int i = 0; i < pData.size() && found != true; i++)
			{
				if(pData.get(i).charAt(j) != ']')
					found = true;
			}

			if(found == false)
			{
				mask[j] = 1;
			}
		}
		return mask;
	}// end align

	/**
	 * This method displays the alignment
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static void displayAlignment(List<Sequence> sData, int numSequences)
	{
		
		List<String> pData = new Vector<String>();

		for (int i = 0; i < Math.min(numSequences+1,sData.size()); i++)
		{
			pData.add(sData.get(i).parseData);
		}
		System.out.println("Displaying alignment ----------------"+pData.size());


		int[] mask = stripDash(pData);

		System.out.println("Alignment from Java parser:");
		for(int j = 0; j < pData.size(); j++)
		{
			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					System.out.print(pData.get(j).charAt(i));
			}		
			if(j == 0)
			{
				System.out.println(" " + sData.get(j).organism);
			}
			else
			{
			System.out.print(" " + sData.get(j).organism + " ");
			for(int x = 0; x < sData.get(j).probablityCount(); x++)
				System.out.print(sData.get(j).getMaxProbability(x) + " ");
			System.out.println();
			}
		}
	}

	/**
	 * This method displays the alignment in FASTA format
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static void displayAlignmentFASTA(List<Sequence> sData, int numSequences)
	{
		
		List<String> pData = new Vector<String>();
		
		for (int i = 0; i < Math.min(numSequences+1,sData.size()); i++)
		{
			pData.add(sData.get(i).parseData);
		}
		System.out.println("Displaying alignment ----------------"+pData.size());

		int[] mask = stripDash(pData);

		System.out.println("Alignment from Java parser:");
		for(int j = 0; j < pData.size(); j++)
		{
			if(j == 0)
			{
				System.out.println("Mask");
			}
			else
			{
				System.out.print(">" + sData.get(j).organism + " ");
				for(int x = 0; x < sData.get(j).probablityCount(); x++)
					System.out.print(sData.get(j).getMaxProbability(x));
				System.out.println();
			}
			for(int i = 0; i < mask.length; i++)
			{
				if (mask[i] == 0)
				{
					char a;
					a = ((String)pData.get(j)).charAt(i);
					System.out.print(a);
				}
			}
			System.out.println();
		}

	}

	/**
	 * This method lists out the alignment of each sequence to the model as a triple
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static void listAlignmentAsCorrespondences(List<Sequence> sData, int numSequences, String SeqName, String ModelName)
	{
		
		List<String> pData = new Vector<String>();
		
		for (int i = 0; i < Math.min(numSequences+1,sData.size()); i++)
		{
			pData.add(sData.get(i).parseData);
		}

		int[] mask = stripDash(pData);

		for(int j = 0; j < pData.size(); j++)
		{
			if(j == 0)
			{
				System.out.println("Mask");
			}
			else
			{
				System.out.print(">" + sData.get(j).organism + " ");
				for(int x = 0; x < sData.get(j).probablityCount(); x++)
					System.out.print(sData.get(j).getMaxProbability(x));
				System.out.println();
			}
			for(int i = 0; i < mask.length; i++)
			{
				if (mask[i] == 0)
				{
					char a;
					a = ((String)pData.get(j)).charAt(i);
					System.out.print(a);
				}
			}
			System.out.println();
		}

	}

	/**
	 * This method displays the alignment in FASTA format
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static void displayMaxLogProbs(List<Sequence> sData, int numSequences, int Motif, int R)
	{		
		System.out.println("Displaying alignment ---------------- "+sData.size());

		System.out.println("Alignment from Java parser:");
		for(int j = 1; j < sData.size(); j++)
		{
			System.out.print(Motif+" "+R+" ");
			for(int x = 0; x < sData.get(j).probablityCount(); x++)
				System.out.print(sData.get(j).getMaxProbability(x) + " ");
			System.out.print(sData.get(j).organism + " ");
			System.out.println();
		}

	}

	
	/**
	 * This method displays the alignment
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static List<String> getAlignment(List<Sequence> sData, int numSequences)
	{
		List<String> pData = new Vector<String>();
		
		for (int i = 0; i < Math.min(sData.size(),numSequences+1); i++)
		{
			pData.add(sData.get(i).parseData);
		}
		
		int[] mask = stripDash(pData);
		String alnm = "";
		Vector<String> alignmentVect = new Vector<String>();

		alnm = "";
		for(int i = 0; i < mask.length; i++)
		{
			if(mask[i] == 0)
				alnm += pData.get(0).charAt(i);
		}		
		alignmentVect.add(alnm);                          // paste in header line
		
		for(int j = 1; j < pData.size(); j++)
		{
			alnm = "";
			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += pData.get(j).charAt(i);
			}		
			alignmentVect.add(alnm);
			alignmentVect.add(">" + sData.get(j).organism);
			alignmentVect.add("Score = " + sData.get(j).getMaxProbability(0));
		}
		return alignmentVect;	
	}
	
	/**
	 * This method displays the alignment
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static List<String> getSortedAlignment(List<String> seqNames, List<String> modNames, int numSequences, int range)
	{
		List<String> alignmentVect = new Vector<String>();
		List<Sequence> sData = new Vector<Sequence>();
		List<String> pData = new Vector<String>();
		double[] modelSums = new double[modNames.size()];
		double[] modelScores = new double[modNames.size()];
		double sum = 0;
		
		// put all sequences in the same Vector
		for(int l = 0; l < seqNames.size(); l++)
		{
			sData.addAll(Alignment.loadFasta(seqNames.get(l)));
		}
		
		// parse all sequences against models
		for(int k = 0; k < modNames.size(); k++)
		{
			sData = Alignment.doParse(sData, numSequences, modNames.get(k), 30);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < sData.get(m).probablityCount(); x++)
			{
				sum = modelSums[x]; // get current sum for this model (x)				
				sum += sData.get(m).getMaxProbability(x);
				modelSums[x] = sum;
			}
		}
		
		System.out.println("Totals for each model: ");
		for(int g = 0; g < modelSums.length; g++)
			System.out.println(modNames.get(g) + " " + modelSums[g]);
		
		System.out.println("Average score for each model unsorted: ");
		for(int g = 0; g < modelSums.length; g++)
		{
			modelScores[g] = (modelSums[g]/sData.size());
			System.out.println(modNames.get(g) + " " + modelScores[g]);
		}
		
		// re-sort models & their totals
		// weird indexing array to keep track of the order things should be in
		// it starts out as [0 1 2 3 4 5]
		// then gets re-sorted according to modelScores
		// ie [0 2 3 1 5 4] means that in other vectors, their index 0 is the first
		// index 3 is second, 2 is the third, ect

		int[] indices = new int[modNames.size()];
		for(int k = 0; k < modNames.size(); k++)
			indices[k] = k;
		
		for(int a = 0; a < modelScores.length; a++)
		   {
	            int max = a; //array position of largest element
	            for(int b = a; b < modelScores.length; b++)
	            {
	                if(modelScores[b] > modelScores[max])
	                	max = b;
	            }
	            // re-sort modelScores or else the sort wouldn't work
	            double dtemp = modelScores[max];
	            modelScores[max] = modelScores[a];
	            modelScores[a] = dtemp;
	            
	            // re-sort indexing array for use with vectors
	            int itemp = indices[max];
	            indices[max] = indices[a];
	            indices[a] = itemp;
		    }
		
		System.out.println("Average score for each model sorted: ");
		for(int g = 0; g < modelSums.length; g++)
		{
			System.out.println(modNames.get(indices[g]) + " " + modelScores[g]);
		}
		
		for (int i = 0; i < sData.size(); i++)
			pData.add(sData.get(i).parseData);

		int[] mask = stripDash(pData);
		String alnm = "";
		
		for(int j = 0; j < pData.size(); j++)
		{
			alnm = "";

			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += pData.get(j).charAt(i);
			}		
			alnm += " " + sData.get(j).organism + " ";
			
			for(int x = 0; x < sData.get(j).probablityCount(); x++)
				alnm += modNames.get(indices[x]) + " score: " + sData.get(j).getMaxProbability(indices[x]) + " ";
			alignmentVect.add(alnm);
		}
		return alignmentVect;
	}

	public static double[] getSortedHLAlignment(List<Sequence> sData, List<String> modNames, int numSequences, int range)
	{
		List<String> alignmentVect = new Vector<String>();
		List<String> pData = new Vector<String>();
		double[] modelSums = new double[modNames.size()];
		double[] modelScores = new double[modNames.size()];
		Vector<String> shortModNames = new Vector<String>();
		Vector<String> tinyModNames = new Vector<String>();
        double[] scores = new double[2*modNames.size()];;

		if(modNames.get(0).contains("http"))
		{
		// remove http:// from model names
			for(int k = 0; k < modNames.size(); k++)
			{
				shortModNames.add(modNames.get(k).substring(33, modNames.get(k).length() - 4));
			}
		}
		else
		{
			shortModNames = new Vector<String>(modNames);
		}

		for (int k=0; k< modNames.size(); k++)
		{
			tinyModNames.add(modNames.get(k).substring(0,8));
		}
		
		// parse all sequences against models
		for(int k = 0; k < modNames.size(); k++)
		{
			sData  = Alignment.doParse(sData, numSequences, modNames.get(k), range);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < sData.get(m).probablityCount(); x++)
			{
				double tempo = sData.get(0).getMaxProbability(x);
				modelSums[x] += tempo;
			}
		}
		
		Formatter fmt = new Formatter(); 

		System.out.println("Average score for each model unsorted: ");
		for(int g = 0; g < modelSums.length; g++)
		{
			modelScores[g] = (modelSums[g]/sData.size());
			System.out.format("%s %12.6f\n",tinyModNames.get(g),modelScores[g]);
			
			fmt = new Formatter();
			scores[g] = Math.max(modelScores[g],-9999);	
		}
		
		// re-sort models & their totals
		// weird indexing array to keep track of the order things should be in
		// it starts out as [0 1 2 3 4 5]
		// then gets re-sorted according to modelScores
		// ie [0 2 3 1 5 4] means that in other vectors, their index 0 is the first
		// index 3 is second, 2 is the third, ect

		int[] indices = new int[modNames.size()];
		for(int k = 0; k < modNames.size(); k++)
			indices[k] = k;
		
		// the following block of code sorts modelScores in an inefficient way
		for(int a = 0; a < modelScores.length; a++)
		   {
	            int max = a; //array position of largest element
	            for(int b = a; b < modelScores.length; b++)
	            {
	                if(modelScores[b] > modelScores[max])
	                	max = b;
	            }
	            // exchange scores
	            double dtemp = modelScores[max];
	            modelScores[max] = modelScores[a];
	            modelScores[a] = dtemp;
	            
	            // exchange indices for use with vectors
	            int itemp = indices[max];
	            indices[max] = indices[a];
	            indices[a] = itemp;
		    }

		System.out.println("Best score for each model, best model first: ");
		for(int g = 0; g < modelSums.length; g++)
		{
			System.out.format("%s %12.6f\n",tinyModNames.get(indices[g]),modelScores[g]);
		}
		
		for (int i = 0; i < sData.size(); i++)
			pData.add(sData.get(i).parseData);

		int[] mask = stripDash(pData);
		String alnm = "";
		
		for(int j = 0; j < pData.size(); j++)
		{
			alnm = "";

			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += pData.get(j).charAt(i);
			}		
			alnm += " " + sData.get(j).organism + " ";
			
			for(int x = 0; x < sData.get(j).probablityCount(); x++)
			{
				fmt = new Formatter();
				fmt.format("%10.6f", sData.get(j).getMaxProbability(indices[x]));
				alnm += tinyModNames.get(indices[x]) + " score: " + fmt + " ";
			}
			alignmentVect.add(alnm);
		}
		return scores;
	}	

	public static double[] getSortedILAlignment(List<Sequence> sData, List<String> modNames, int numSequences, int range)
	{
		Vector<String> pData = new Vector<String>();                           // parse data
		double[] modelSums = new double[modNames.size()];      // sum of alignment scores
		double[] rmodelSums = new double[modNames.size()];     // sum with sequences reversed
		double[] modelScores = new double[modNames.size()];
		double[] rmodelScores = new double[modNames.size()];
		List<String> shortModNames = new Vector<String>();                   // for easier display
		List<String> tinyModNames = new Vector<String>();                    // for even easier display
		int[] reversed = new int[modNames.size()];             // is best model reversed?
 //       String scores = "";
        double[] scores = new double[2*modNames.size()];       // all scores computed

        List<Sequence> rsData = Alignment.reverse(numSequences, sData);  // reversed sequence data
		
		if(((String)modNames.get(0)).contains("http"))         // look online for models
		{
		// remove http:// from model names
			for(int k = 0; k < modNames.size(); k++)
			{
				shortModNames.add(modNames.get(k).substring(33, modNames.get(k).length()-4));
			}
		}
		else
		{
			shortModNames = new Vector<String>(modNames);
		}

		for (int k=0; k< modNames.size(); k++)
		{
			tinyModNames.add(modNames.get(k).substring(0,6));
		}
		
		// parse sequence data in sData against models
		for(int k = 0; k < modNames.size(); k++)
		{
			sData  = Alignment.doParse(sData, numSequences, modNames.get(k), range);
			rsData = Alignment.doParse(rsData, numSequences, modNames.get(k), range);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < sData.get(m).probablityCount(); x++)
			{
				// there is no good reason for having to do these crazy manipulations
				// in order to get the value of a double variable, but at least this works
			    double tempo = sData.get(m).getMaxProbability(x);   
				modelSums[x] += tempo;

				tempo = rsData.get(m).getMaxProbability(x);
				rmodelSums[x] += tempo;
			}
		}
		
		Formatter fmt = new Formatter(); 
		
		for(int g = 0; g < modelSums.length; g++)
		{
			modelScores[g] = (modelSums[g]/(sData.size()-1));
			rmodelScores[g] = (rmodelSums[g]/(sData.size()-1));
			
			fmt = new Formatter();
			scores[2*g]   = Math.max(modelScores[g],-9999);
			scores[2*g+1] = Math.max(rmodelScores[g],-9999);
			
			if(rmodelScores[g] > modelScores[g])          // choose between forward and reversed for each model
			{
				modelScores[g] = rmodelScores[g];
				shortModNames.set(g, shortModNames.get(g)+" reversed");
				tinyModNames.set(g, tinyModNames.get(g)+" reversed");
                reversed[g] = 1;
			}
			else {
				reversed[g] = 0;
			}
		}
		
		// re-sort models & their totals
		// weird indexing array to keep track of the order things should be in
		// it starts out as [0 1 2 3 4 5]
		// then gets re-sorted according to modelScores
		// ie [0 2 3 1 5 4] means that in other vectors, their index 0 is the first
		// index 3 is second, 2 is the third, ect

		int[] indices = new int[modNames.size()];
		for(int k = 0; k < modNames.size(); k++)
			indices[k] = k;
		
		// the following block of code sorts modelScores in an inefficient way
		for(int a = 0; a < modelScores.length; a++)
		   {
	            int max = a; //array position of largest element
	            for(int b = a; b < modelScores.length; b++)
	            {
	                if(modelScores[b] > modelScores[max])
	                	max = b;
	            }
	            // exchange scores
	            double dtemp = modelScores[max];
	            modelScores[max] = modelScores[a];
	            modelScores[a] = dtemp;
	            
	            // exchange indices for use with vectors
	            int itemp = indices[max];
	            indices[max] = indices[a];
	            indices[a] = itemp;
		    }
		
		System.out.println("");
		System.out.println("Best score for each model, best model first: ");
		for(int g = 0; g < modelSums.length; g++)
		{
			if (reversed[indices[g]] == 0)
				System.out.format("%s          %12.6f\n",tinyModNames.get(indices[g]),modelScores[g]);
			else
				System.out.format("%s %12.6f\n",tinyModNames.get(indices[g]),modelScores[g]);
		}
		
		//		 align again to get max prob alignment
        if (reversed[indices[0]]==0) {
        	sData  = Alignment.doParse(sData, numSequences, modNames.get(indices[0]), range);
        	System.out.println("parsing forward again " + modNames.get(indices[0]) + " " + indices[0]);
    		for (int i = 0; i < sData.size(); i++)
    			pData.add(sData.get(i).parseData);  
        }
        else {
        	rsData = Alignment.doParse(rsData, numSequences, modNames.get(indices[0]),range);
        	System.out.println("parsing reversed again " + modNames.get(indices[0]) + " " + indices[0]);
    		for (int i = 0; i < sData.size(); i++)
    			pData.add(rsData.get(i).parseData); 
        }
		
		

		int[] mask = stripDash(pData);
		String alnm = "";

		for(int j = 0; j < Math.min(pData.size(),30); j++)    // at most 30 sequences from each
		{
			alnm = "";

			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += pData.get(j).charAt(i);
			}		
			
			if (j > 0) {
			alnm += " "+((Sequence)sData.get(j)).organism+" ";
			for(int x = 0; x < modNames.size(); x++)
			{
				fmt = new Formatter();
				if (reversed[indices[x]] == 0)
					fmt.format("%12.6f", sData.get(j).getMaxProbability(indices[x]));
				else
					fmt.format("%12.6f", rsData.get(j).getMaxProbability(indices[x]));
				alnm += tinyModNames.get(indices[x]) + " score: " + fmt + " ";
			}
			}
        	System.out.println(alnm);
		}

		return scores;
	}	

	public static void printAlignment(List<String> aData, int numChars)
	{
		for(int f = 0; f < aData.get(1).length()/numChars+1; f++)
		{
			for(int j = 0; j < aData.size(); j++)
			{
				if((f+1)*numChars > aData.get(j).length())
				{
					System.out.println(aData.get(j).substring(f*numChars, aData.get(j).length()));
				}
				else
					System.out.println(aData.get(j).substring(f*numChars, (f+1)*numChars));
			}
			System.out.println();
		}
	}

	public static String tempdisplayAlignmentFASTA(List<Sequence> sData, int numSequences)
	{
		String temp = "";
		Vector<String> pData = new Vector<String>();
		
		for (int i = 0; i < Math.min(numSequences+1,sData.size()); i++)
		{
			pData.add(sData.get(i).parseData);
		}
		temp += "Displaying alignment ----------------"+pData.size()+"\n";


		int[] mask = stripDash(pData);
		System.out.println("Alignment from Java parser:");
		for(int j = 0; j < pData.size(); j++)
		{
			if(j == 0)
			{
				temp += "Mask\n";
			}
			else
			{
				temp += ">" + sData.get(j).organism + " ";
				for(int x = 0; x < sData.get(j).probablityCount(); x++)
					temp += sData.get(j).getMaxProbability(x);
				temp += "\n";
			}
			for(int i = 0; i < mask.length; i++)
			{
				if (mask[i] == 0)
				{
					char a;
					a = ((String)pData.get(j)).charAt(i);
					if ((a!='{') && (a!='}') && (a!='<') && (a!='>')) 
						temp += a;
				}
			}
			temp += "\n";
		}

		return temp;
	}
	
	public static ParseData doParse2(List<Sequence> sData, int numSequences, String nodeFileName, int range)
	{
		Node current;
		List<Double> mProbs = new Vector<Double>();
		List<List<Double>> probsM = new Vector<List<Double>>();
		
		Sequence S = new Sequence("","");                        // blank sequence to use repeatedly 
		S.addNodeData(nodeFileName);	                         // only add node data once

		sData.get(0).parseData = ((InitialNode)S.first).header();

		Sequence firstS = sData.get(1);          // first sequence, which matches the model
		firstS.setNucleotides();
		firstS.setArrays();
		
		for (int i = 1; i < Math.min(sData.size(),numSequences+1); i++) 
		{
			S.organism = sData.get(i).organism;             // focus on one sequence
			S.letters  = sData.get(i).letters;
			
			// System.out.println("Alignment.doParse: " + S.letters);
			
			S.setNucleotides();                                    // strip dashes from sequence
			S.setArrays();                                         // define cti, itc, convert letters to numbers

			// keep copies of cti and itc for the first sequence, which must match the model
			S.ctiFirst = new int[firstS.cti.length];            
			S.itcFirst = new int[firstS.itc.length];
			
			for (int j = 0; j < S.ctiFirst.length; j++)
				S.ctiFirst[j] = firstS.cti[j];
			for (int j = 0; j < S.itcFirst.length; j++)
				S.itcFirst[j] = firstS.itc[j];

			S.parseSequence(range);                             	      // parse this sequence

			// We also need to store the parse information somewhere!  All we have is a parse sequence.
			// We need to store the max log probabilities too

			mProbs = new Vector<Double>();
			current = S.first;
			while(current != null)
			{
				mProbs.add(new Double(current.optimalMaxLogProb));
				current =  current.next;
			}
			sData.get(i).appendProbabilities(mProbs);
			if (S.first.optimalMaxLogProb < -99999999)
			{
				current = S.first.next;
				while (current != null)
				{
					if (current.previous.currentMaxLogProb < -99999999 && current.currentMaxLogProb > -999999999)
					{
					}
					current = current.next;				
				} // TODO REMOVE
			}
			probsM.add(mProbs);
			sData.get(i).parseData = ((InitialNode)S.first).showParse(S.nucleotides);
			
		}
		Double[][] probsArray = Alignment.vec2array(probsM);
		ParseData PD = new ParseData();
		PD.sdata = sData;
		PD.probsM = probsArray;
		return PD;
	}

	public static Double[][] vec2array(List<List<Double>> probsM)
	{
		int n = probsM.size();
		List<Double> dummy = probsM.get(0);
		int m = dummy.size();
		Double[][] probsArray = new Double[n][m];
		for(int i=0; i<n; i++)
		{
			List<Double> vec = probsM.get(i);
			for(int j=0; j<m; j++)
			{	
				probsArray[i][j] = vec.get(j);
			}
		}
		return probsArray;
	}

	public static double[][] getILScores(List<Sequence> sData, List<String> modNames, int numSequences, int range)
	{
        double[][] scores = new double[sData.size()][modNames.size()];       // all scores computed
		
		// parse sequence data in sData against models
		for(int k = 0; k < modNames.size(); k++)
		{
			sData  = Alignment.doParse(sData, numSequences, modNames.get(k), range);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < sData.get(m).probablityCount(); x++)
			{
				// there is no good reason for having to do these crazy manipulations
				// in order to get the value of a double variable, but at least this works
//				String temp = String.valueOf(((Vector)sData.elementAt(m).maxLogProbs.get(x)).get(0)); // get score for this sequence(m)
				double tempo = sData.get(m).getMaxLogProbability(x, 0);   
				scores[m][x] = tempo;
			}
		}
		return scores;
	}
	public static double[] getILScoresSingle(List<Sequence> sData, String modName, int numSequences, int range)
	{
		double[] scores = new double[sData.size()-1];       // all scores computed
		sData  = Alignment.doParse(sData,numSequences,modName,range);
		// add up model scores for each sequence
		for(int m = 0; m < sData.size()-1; m++)
		{
				double tempo = sData.get(m+1).getMaxLogProbability(0, 0);   
				scores[m] = tempo;
		}
		return scores;
	}
	
	//Takes a JAR3D query and submits results to MySQL database
	public static List<LoopResult> doILdbQuery(int loopID, Query query, List<Sequence> sData, List<String> modNames, HashMap<String, MotifGroup> groupData, int numSequences, int range)
	{
		double[] modelSums = new double[modNames.size()];      // sum of alignment scores
		double[] rmodelSums = new double[modNames.size()];     // sum with sequences reversed
		double[] modelScores = new double[modNames.size()];
		double[] rmodelScores = new double[modNames.size()];
		double[][] modelScoreMat = new double[modNames.size()][sData.size()];
		double[][] rmodelScoreMat = new double[modNames.size()][sData.size()];
		List<String> shortModNames = new Vector<String>();                   // for easier display
		int[] reversed = new int[modNames.size()];             // is best model reversed?
		double[] scores = new double[2*modNames.size()];
		List<LoopResult> loopRes = new Vector<LoopResult>();
		
		List<Sequence> rsData = Alignment.reverse(numSequences, sData);  // reversed sequence data
	    
	    shortModNames = new Vector<String>(modNames);
	    Vector<String> tinyModNames = new Vector<String>(shortModNames);
//Parse all sequences against all groups	    
	    MotifGroup group;
	    for(int k = 0; k < modNames.size(); k++)
		{
	    	group = groupData.get(modNames.get(k));
			sData  = Alignment.doParse(sData,numSequences,group.Model,range,Boolean.TRUE);
			rsData = Alignment.doParse(rsData,numSequences,group.Model,range,Boolean.TRUE);
		}
		
//Add up model scores for each sequence, find mean score, compare regular and reversed scores
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < sData.get(m).probablityCount(); x++)
			{
				double tempo = sData.get(m).getMaxLogProbability(x, 0);   
				modelSums[x] += tempo;
				modelScoreMat[x][m] = tempo;

				tempo = sData.get(m).getMaxLogProbability(x, 0);
				rmodelSums[x] += tempo;
				rmodelScoreMat[x][m] = tempo;
			}
		}
		for(int g = 0; g < modelSums.length; g++)
		{
			modelScores[g] = (modelSums[g]/(sData.size()-1));
			rmodelScores[g] = (rmodelSums[g]/(sData.size()-1));
			scores[2*g]   = Math.max(modelScores[g],-9999);
			scores[2*g+1] = Math.max(rmodelScores[g],-9999);
		
			if(rmodelScores[g] > modelScores[g])          // choose between forward and reversed for each model
			{
				modelScores[g] = rmodelScores[g];
				shortModNames.set(g, shortModNames.get(g)+" reversed");
		        reversed[g] = 1;
			}
			else {
				reversed[g] = 0;
			}
		}
		
		// re-sort models & their totals
		// weird indexing array to keep track of the order things should be in
		// it starts out as [0 1 2 3 4 5]
		// then gets re-sorted according to modelScores
		// ie [0 2 3 1 5 4] means that in other vectors, their index 0 is the first
		// index 3 is second, 2 is the third, ect

		int[] indices = new int[modNames.size()];
		for(int k = 0; k < modNames.size(); k++)
			indices[k] = k;
		

//Calculate extra information (quantiles, edit distances) and out put
		int numInputSeqs = sData.size()-1;
		for(int g = 0; g < modNames.size(); g++)
		{
			int index = indices[g];
			String groupName = tinyModNames.get(index);
			String sig;
			boolean rev;
			double[] groupScores = new double[sData.size()];
			group = groupData.get(groupName);
			if (reversed[index] == 0){  //not reversed
				rev = Boolean.FALSE;
				sig = group.Signature[0]; 
				for(int col = 1; col < numInputSeqs+1; col++)
				{
					groupScores[col-1] = modelScoreMat[index][col];
				}
			}
			else{  //reversed
				rev = Boolean.TRUE;
				sig = group.Signature[1];
				for(int col = 1; col < numInputSeqs+1; col++)
				{
					groupScores[col-1] = rmodelScoreMat[index][col];
				}
			}
			//Calculate quantiles
			double[] quants = webJAR3D.getQuantilesA(groupScores, group);
			//Calculate edit distances
			Vector<Sequence> modsData = Alignment.parseFastaText(group.Sequences,0,0);
			int[][] EditDistances = SimpleAlign.calcILEditDistances(sData,modsData,rev);
			int[] minDist = new int[EditDistances.length];
			for(int i =0; i < EditDistances.length; i ++){
				minDist[i] = ArrayMath.min(EditDistances[i]);
			}
			if(Boolean.TRUE){   //put goodness of fit checks here later?
				Vector<SequenceResult> seqRes = new Vector<SequenceResult>();
				for(int m = 0; m < sData.size() - 1; m++)
				{
					MutableSequenceResults seqR = new MutableSequenceResults(groupName,groupScores[m],quants[m],
							minDist[m],rev,loopID);
					seqR.setQuery(query);
					seqR.setSequenceId(String.valueOf(m));
					seqRes.add(seqR);
				}
				ImmutableLoopResult loopR = new ImmutableLoopResult(loopID, groupName, rev, sig, seqRes, "NA", "");
				loopR.setQuery(query);
				loopRes.add(loopR);
			}
		}
		return loopRes;
	}
	
	//Overloaded doParse that can take group data as a string instead of a file name
	public static List<Sequence> doParse(List<Sequence> sData, int numSequences, String nodeInfo, int range, boolean fullModelText)
	{
		Node current;
		List<Double> mProbs = new Vector<Double>();

		Sequence S = new Sequence("","");                    // blank sequence to use repeatedly
		if(fullModelText){
			S.addNodeDataModelText(nodeInfo);				 // add model data from string
		}else{
			S.addNodeData(nodeInfo);	                     // add model data from file
		}

		sData.get(0).parseData = ((InitialNode)S.first).header(); // add a header line

		Sequence firstS = sData.get(1);          // first sequence, which matches the model
		firstS.setNucleotides();
		firstS.setArrays();
		
		for (int i = 1; i < Math.min(sData.size(),numSequences+1); i++) 
		{
			S.organism = sData.get(i).organism;             // focus on one sequence
			S.letters  = sData.get(i).letters;			
			S.setNucleotides();                                    // strip dashes from sequence
			S.setArrays();                                         // define cti, itc, convert letters to numbers

			// keep copies of cti and itc for the first sequence, which must match the model
			S.ctiFirst = new int[firstS.cti.length];            
			S.itcFirst = new int[firstS.itc.length];
			
			for (int j = 0; j < S.ctiFirst.length; j++)
				S.ctiFirst[j] = firstS.cti[j];
			for (int j = 0; j < S.itcFirst.length; j++)
				S.itcFirst[j] = firstS.itc[j];

			S.parseSequence(range);                             	      // parse this sequence

			// We also need to store the parse information somewhere!  All we have is a parse sequence.
			// We need to store the max log probabilities too

			mProbs = new Vector<Double>();
			current = S.first;
			while(current != null)
			{
				mProbs.add(new Double(current.optimalMaxLogProb));
				current =  current.next;
			}
			sData.get(i).appendProbabilities(mProbs);
			if (S.first.optimalMaxLogProb < -99999999)
			{
				current = S.first.next;
				while (current != null)
				{
					if (current.previous.currentMaxLogProb < -99999999 && current.currentMaxLogProb > -999999999)
					{
					}
					current = current.next;				
				}
			}
			
			sData.get(i).parseData = ((InitialNode)S.first).showParse(S.nucleotides);
			
			String correspondences = ((InitialNode)S.first).showCorrespondences(S.nucleotides);

			correspondences = correspondences.replace("JAR3D_aligns_to", "aligns_to_JAR3D");
			
			String SF = "Sequence_"+i+"_"+S.organism;
			SF = SF.replace(" ","_");
			correspondences = correspondences.replace("SSS",SF);
			
			String NF = "MMM";                           // model name goes here
			correspondences = correspondences.replace("MMM", NF);
			sData.get(i).correspondences = correspondences;
		}
		return sData;
	}

}
