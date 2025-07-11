package edu.bgsu.rna.jar3d;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Formatter;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.results.BasicLoopResult;
import edu.bgsu.rna.jar3d.results.LoopResult;
import edu.bgsu.rna.jar3d.results.BasicSequenceResult;
import edu.bgsu.rna.jar3d.results.SequenceResult;

/**
 * This is the Alignment class, it has methods that are used to facilitate aligning sequences
 * @author meg pirrung
 *
 */
public class Alignment {

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
	public static Vector<Sequence> loadFastaColumnsDNA(String fileName, int StartCol, int EndCol, int DNA) 		{
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
				// Please retain the following commented-out line for debugging
				// System.out.println("Reading sequence file from relative path");
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
					if (DNA > 0) {
						// change A to U, G to C, etc.
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

			// System.out.println(String.format("Alignment.loadFastaColumnsDNA: filename %s organism %s letters %s",fileName,organism,letters));

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
	 * @param sequenceData
	 * @return
	 */
	public static List<Sequence> reverse(List<Sequence> sequenceData)
	{
		List<Sequence> reverseSData = new Vector<Sequence>();

		String temp="";
		String organism="";
		String letters="";
		String second = "";
		String first  = "";

		reverseSData.add(new Sequence("","")); // save space for header information to be added later as element 0 of this vector

		for (int i=1; i < sequenceData.size(); i++)
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
	 * This method rotates strand orders r times in sequenceData and returns a new vector of sequence data
	 * @param sequenceData
	 * @return
	 */
	public static List<Sequence> rotate(List<Sequence> sequenceData, int rotation) {
		List<Sequence> rotatedSData = new Vector<Sequence>();

		String temp="";
		String organism="";
		String letters= "";
		String second = "";
		String first  = "";

		rotatedSData.add(new Sequence("","")); // save space for header information to be added later as element 0 of this vector

		for (int s=1; s < sequenceData.size(); s++)		{
			temp = sequenceData.get(s).letters;
			organism = sequenceData.get(s).organism;

			String[] parts = temp.split("\\*");  // split on "*"
			int n = parts.length;
			StringBuilder result = new StringBuilder();

			for (int k = 0; k < n; k++) {
				result.append(parts[(k + rotation) % n]);
				if (k < n-1) {
					result.append("*");
				}
			}

		    letters = result.toString();

			System.out.println(String.format("Alignment rotate: rotation %d",rotation));
			System.out.println(String.format("Alignment rotate: organism %s",organism));
			System.out.println(String.format("Alignment rotate: letters  %s",letters));

			// print the stack so we know how we got here
			// StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
			// StringBuilder sb = new StringBuilder();
			// for (StackTraceElement element : stackTrace) {
			// 	System.out.println(element.toString());
			// }

			// pause so we can see the stack trace
			// try {
			// 	Thread.sleep(1000);
			// } catch (InterruptedException e) {
			// 	e.printStackTrace();  // or handle interruption another way
			// }

			rotatedSData.add(new Sequence(organism,letters));

		}
		return rotatedSData;
	}


	/**
	 * This method initiates parsing.
	 * @param sData
	 * @param nodeFileName
	 * @param range
	 * @return
	 */
	public static List<Sequence> doParse(List<Sequence> sData, String nodeFileName, int range)
	{
		return Alignment.doParse(sData, nodeFileName, range, false, false);
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
	 * This method displays the alignment in FASTA format
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static void displayAlignmentFASTA(List<Sequence> sData)
	{

		List<String> pData = new Vector<String>();

		for (int i = 0; i < sData.size(); i++)
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
				for(int x = 0; x < sData.get(j).getMaxNodeLogProbabilitySize(); x++)
					System.out.print(sData.get(j).getMaxNodeLogProbability(x));
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
	public static void listAlignmentAsCorrespondences(List<Sequence> sData, String SeqName, String ModelName)
	{

		List<String> pData = new Vector<String>();

		for (int i = 0; i < sData.size(); i++)
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
				for(int x = 0; x < sData.get(j).getMaxNodeLogProbabilitySize(); x++)
					System.out.print(sData.get(j).getMaxNodeLogProbability(x));
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

	public static double[] getSortedHLAlignment(List<Sequence> sData, List<String> modNames, int range)
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
			sData  = Alignment.doParse(sData, modNames.get(k), range);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < sData.get(m).getMaxNodeLogProbabilitySize(); x++)
			{
				double tempo = sData.get(0).getMaxNodeLogProbability(x);
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

			for(int x = 0; x < sData.get(j).getMaxNodeLogProbabilitySize(); x++)
			{
				fmt = new Formatter();
				fmt.format("%10.6f", sData.get(j).getMaxNodeLogProbability(indices[x]));
				alnm += tinyModNames.get(indices[x]) + " score: " + fmt + " ";
			}
			alignmentVect.add(alnm);
		}

		fmt.close();

		return scores;
	}

	public static double[] getSortedILAlignment(List<Sequence> sData, List<String> modNames, int range) 	{
		Vector<String> pData = new Vector<String>();                           // parse data
		double[] modelSums = new double[modNames.size()];      // sum of alignment scores
		double[] rmodelSums = new double[modNames.size()];     // sum with sequences reversed
		double[] modelScores = new double[modNames.size()];
		double[] rmodelScores = new double[modNames.size()];
		List<String> shortModNames = new Vector<String>();                   // for easier display
		List<String> tinyModNames = new Vector<String>();                    // for even easier display
		int[] reversed = new int[modNames.size()];             // is best model reversed?
        double[] scores = new double[2*modNames.size()];       // all scores computed

        List<Sequence> rsData = Alignment.reverse(sData);  // reversed sequence data

		if(modNames.get(0).contains("http"))         // look online for models
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
			sData  = Alignment.doParse(sData, modNames.get(k), range);
			rsData = Alignment.doParse(rsData, modNames.get(k), range);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < sData.get(m).getMaxNodeLogProbabilitySize(); x++)
			{
				// there is no good reason for having to do these crazy manipulations
				// in order to get the value of a double variable, but at least this works
			    double tempo = sData.get(m).getMaxNodeLogProbability(x);
				modelSums[x] += tempo;

				tempo = rsData.get(m).getMaxNodeLogProbability(x);
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
        	sData  = Alignment.doParse(sData, modNames.get(indices[0]), range);
        	System.out.println("parsing forward again " + modNames.get(indices[0]) + " " + indices[0]);
    		for (int i = 0; i < sData.size(); i++)
    			pData.add(sData.get(i).parseData);
        }
        else {
        	rsData = Alignment.doParse(rsData, modNames.get(indices[0]),range);
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
					fmt.format("%12.6f", sData.get(j).getMaxNodeLogProbability(indices[x]));
				else
					fmt.format("%12.6f", rsData.get(j).getMaxNodeLogProbability(indices[x]));
				alnm += tinyModNames.get(indices[x]) + " score: " + fmt + " ";
			}
			}
        	System.out.println(alnm);
		}

		fmt.close();

		return scores;
	}


	// doParse2 is probably outmoded and should be avoided - CLZ 2013-03-04

	public static ParseData doParse2(List<Sequence> sData, String nodeFileName, int range)
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

		for (int i = 1; i < sData.size(); i++)
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
			sData.get(i).appendNodeProbabilities(mProbs);
			probsM.add(mProbs);
			sData.get(i).parseData = ((InitialNode)S.first).showParse(S.nucleotides);

		}
		Double[][] probsArray = Alignment.vec2array(probsM);
		ParseData PD = new ParseData();
		PD.sdata = sData;
		PD.probsM = probsArray;
		return PD;
	}

	public static Double[][] vec2array(List<List<Double>> probsM) {
		int n = probsM.size();
		Double[][] probsArray = new Double[n][];
        Double[] blank = new Double[0];
		for(int i=0; i<n; i++) {
            probsArray[i] = probsM.get(i).toArray(blank);
		}
		return probsArray;
	}

	public static double[][] getILScores(List<Sequence> sData, List<String> modNames, int range)
	{
        double[][] scores = new double[sData.size()][modNames.size()];       // all scores computed

		// parse sequence data in sData against models
		for(int k = 0; k < modNames.size(); k++)
		{
			sData  = Alignment.doParse(sData, modNames.get(k), range);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < sData.get(m).getMaxNodeLogProbabilitySize(); x++)
			{
				// there is no good reason for having to do these crazy manipulations
				// in order to get the value of a double variable, but at least this works
				double tempo = sData.get(m).getMaxNodeLogProbabilityOf(x, 0);
				scores[m][x] = tempo;
			}
		}
		return scores;
	}

	public static double[] getILScoresSingle(List<Sequence> sData, String modName, int range) {
		// score sequences in sData against one model
		// not just for IL

		double[] scores = new double[sData.size()-1];       // all scores computed, starting at 0
		sData  = Alignment.doParse(sData, modName, range);  // entry 0 for header, start at 1

		// add up model scores for each sequence
		for (int s = 1; s < sData.size(); s++) {
			// double tempo = sData.get(s).getMaxNodeLogProbabilityOf(0, 0);  // stopped working! same number for all s, even though fine in doParse
			double tempo = sData.get(s).getMaxLogProbability(0);
			scores[s-1] = tempo;

			// System.out.println("getILScoresSingle:  getMaxNodeLogProbabilityOf: "+sData.get(s).getMaxNodeLogProbabilityOf(0, 0));
			// System.out.println("getILScoresSingle:  scores[m]:                  "+scores[s-1]);
			// System.out.println("getILScoresSingle:  getMaxLogProbability:       "+sData.get(s).getMaxLogProbability(0));
		}
		return scores;
	}
	public static List<LoopResult> doLoopDBQuery(Loop loop, List<String> modNames, HashMap<String, MotifGroup> groupData,
			int range) {
		return doLoopDBQuery(loop, modNames, groupData, range, false);
	}

	public static List<LoopResult> doLoopDBQuery(Loop loop, List<String> modNames, HashMap<String, MotifGroup> groupData,
			int range, boolean saveSeqRes) {

		List<Sequence> sData = loop.getSequences();

		Query query = loop.getQuery();

		// 2013-11-05 CLZ Last argument cannot be null, must be something like "IL"
		// 2013-11-05 CLZ Last argument seems to be reliably obtained from loop
		List<LoopResult> results = doLoopDBQuery((int)loop.getId(), query, sData, modNames, groupData, range, loop, loop.getTypeString(), saveSeqRes);

		for(LoopResult result: results) {
			result.setLoop(loop);
		}

		return results;
	}

	public static List<SequenceResult> doSingleDBQuery(Loop loop, MotifGroup group, String model_id, int range) {

		List<Sequence> sData = loop.getSequences();
		Query query = loop.getQuery();

		// 2013-11-05 CLZ Last argument cannot be null, must be something like "IL"
		// 2013-11-05 CLZ Last argument seems to be reliably obtained from loop
		List<SequenceResult> results = doSingleDBQuery(query, sData, group, model_id, range, loop, loop.getTypeString());

		return results;
	}

	//Takes a JAR3D query and submits results to MySQL database
		public static List<LoopResult> doLoopDBQuery(int loopID, Query query, List<Sequence> sData, List<String> modNames,
				HashMap<String, MotifGroup> groupData, int range, Loop loop, String type) {
			return doLoopDBQuery(loopID, query, sData, modNames,
					groupData, range, loop, type, false);
		}

	//Takes a JAR3D query and submits results to MySQL database
	public static List<LoopResult> doLoopDBQuery(int loopID, Query query, List<Sequence> sData, List<String> modNames,
			HashMap<String, MotifGroup> groupData, int range, Loop loop, String type, Boolean saveSeqRes) {

		// new code to handle all required rotations at once
		int maxRotation = 0;
		if (type.equalsIgnoreCase("HL")) {
			maxRotation = 0;
		} else if (type.equalsIgnoreCase("IL")){
			maxRotation = 1;
		} else if (type.startsWith("J")) {
			// J3, J4, ...
			maxRotation = Integer.parseInt(type.substring(1)) - 1;
		}

		List<Sequence>[] rotatedData = new Vector[maxRotation+1];  // Array of Vector<Sequence>
		int[] bestRotation = new int[modNames.size()];             // 0 for HL, 0 or 1 for IL, etc.
		double[] scores = new double[modNames.size()];
		List<LoopResult> loopRes = new Vector<LoopResult>();

		// rotate sequence data however many times are needed; 0 for HL, 0 and 1 for IL, etc.
		for (int r = 0; r <= maxRotation; r++) {
			rotatedData[r] = Alignment.rotate(sData,r);
		}

		List<String> shortModNames = new Vector<String>();                   // for easier display
	    shortModNames = new Vector<String>(modNames);
		Vector<String> tinyModNames = new Vector<String>(shortModNames);

		// Parse all sequences against all groups in all rotations
	    MotifGroup group;
	    for (int k = 0; k < modNames.size(); k++) {
	    	group = groupData.get(modNames.get(k));
			for (int r = 0; r <= maxRotation; r++) {
				rotatedData[r] = Alignment.doParse(rotatedData[r], group.Model, range, true, false, false);
			}
		}

	    // Record model scores for each sequence, sum over sequences
		double[][] modelRotationScoreSum = new double[modNames.size()][maxRotation+1];      // sum of alignment scores
		double[][][] modelRotationSequenceScore = new double[modNames.size()][maxRotation+1][sData.size()];

		for (int r = 0; r <= maxRotation; r++) {
			for (int m = 0; m < modNames.size(); m++) {
				for (int s = 1; s < sData.size(); s++) {
					double tempo = rotatedData[r].get(s).getMaxLogProbability(m);
					modelRotationScoreSum[m][r] += tempo;
					modelRotationSequenceScore[m][r][s] = tempo;
				}
			}
		}

		// for each group find the best rotation
	    for (int m = 0; m < modNames.size(); m++) {
			scores[m] = -9999;
			bestRotation[m] = 0;
			for (int r = 0; r <= maxRotation; r++) {
				// number of sequences is sData.size()-1, 0th place is empty to leave space for headers
				double tempScore = modelRotationScoreSum[m][r]/(sData.size()-1);
				if (tempScore > scores[m]) {
					scores[m] = tempScore;
					bestRotation[m] = r;
				}
			}
			System.out.println(String.format("Alignment: groupName %s best rotation %d score %8.2f",modNames.get(m),bestRotation[m],scores[m]));
		}

		// re-sort models & their totals ... doesn't seem to be used anymore
		// weird indexing array to keep track of the order things should be in
		// it starts out as [0 1 2 3 4 5]
		// then gets re-sorted according to modelScores
		// ie [0 2 3 1 5 4] means that in other vectors, their index 0 is the first
		// index 3 is second, 2 is the third, ect
		int[] indices = new int[modNames.size()];
		for(int m = 0; m < modNames.size(); m++) {
			indices[m] = m;
		}

		// Calculate extra information (edit distances, cutoffs) and output
		for (int m = 0; m < modNames.size(); m++) {

			int index = indices[m];
			String groupName = tinyModNames.get(index);
			String sig;
			int rotation = bestRotation[m];
			double[] groupScores = new double[sData.size()];
			group = groupData.get(groupName);

			// System.out.println(String.format("Alignment: groupName %s",groupName));

			sig = group.Signature[0]; // don't reverse the signature, too confusing

			System.out.println(String.format("group.Signature %s", group.Signature));

			for (int s = 1; s < sData.size(); s++) {
				groupScores[s-1] = modelRotationSequenceScore[index][rotation][s];
			}

			double[] deficits = getDeficits(group.Best_Score, groupScores);

			// sequences are already rotated; pass in the optimal rotation, 0 additional rotation
			// EditDistance functions need to work on J3, J4, etc.
			int[] InteriorMinDist = getMinEditDistance(group,rotatedData[rotation],type,0,true);
			int[] FullMinDist = getMinEditDistance(group,rotatedData[rotation],type,0,false);

			boolean[] cutoffs = getCutoffs(group,groupScores,InteriorMinDist);
			double[] cutoffscores = getCutoffScores(group,groupScores,InteriorMinDist);

			if (true) {
				LoopResult loopR;
				if (saveSeqRes) {
					List<SequenceResult> seqRes = new ArrayList<SequenceResult>();
					for (int s = 1; s < sData.size(); s++) {
						SequenceResult seqR = new BasicSequenceResult(sData.get(s), groupScores[s-1], deficits[s-1], InteriorMinDist[s-1], FullMinDist[s-1],rotation,cutoffs[s-1],cutoffscores[s-1]);
						seqRes.add(seqR);
					}
					loopR = new BasicLoopResult(groupName, rotation, sig, seqRes);
				} else {
					double medianScore = ArrayMath.median(groupScores);
					double meanScore = ArrayMath.mean(groupScores);
					double meanDeficit = ArrayMath.mean(deficits);
					double medianDeficit = ArrayMath.median(deficits);
					double meanInteriorEditDistance = ArrayMath.mean(InteriorMinDist);
					double medianInteriorEditDistance = ArrayMath.median(InteriorMinDist);
					double meanFullEditDistance = ArrayMath.mean(FullMinDist);
					double medianFullEditDistance = ArrayMath.median(FullMinDist);
					double meanCutoff = ArrayMath.mean(cutoffs);
					double meanCutoffScore = ArrayMath.mean(cutoffscores);

					loopR = new BasicLoopResult(groupName, rotation, sig,
						medianScore, meanScore,
						meanDeficit, medianDeficit,
						meanInteriorEditDistance, medianInteriorEditDistance,
						meanFullEditDistance, medianFullEditDistance,
						meanCutoff, meanCutoffScore);
				}
				loopR.setLoop(loop);
				loopRes.add(loopR);
			}
		}
		//Try to clean up the memory we used
		System.gc();
		System.runFinalization();
		System.gc();

		return loopRes;
	}

	// Takes a JAR3D query and submits results to MySQL database
	public static List<SequenceResult> doSingleDBQuery(Query query, List<Sequence> sData,
			MotifGroup group, String model_id, int range, Loop loop, String type) {

		// TODO 2013-11-07 CLZ generalize this so that it can apply to HL, IL, 3WJ, 4WJ, etc.
		// Currently there is only space for rotations 0 and 1
		double modelSum=0;
		double rmodelSum=0;
		double modelScore=-9999;
		double rmodelScore=-9999;
		double[] modelScores = new double[sData.size()];
		double[] rmodelScores = new double[sData.size()];
		int reversed;
		List<Sequence> rsData = new Vector<Sequence>();

		if(type.equalsIgnoreCase("IL")){
			rsData = Alignment.reverse(sData);  // reversed sequence data
		}

	    sData  = Alignment.doParse(sData, group.Model, range, true, true, false);
			if(type.equalsIgnoreCase("IL")){
				rsData = Alignment.doParse(rsData, group.Model, range, true, true, false);
			}
	    //Add up model scores for each sequence, find mean score, compare regular and reversed scores
	    for(int m = 1; m < sData.size(); m++)
	    {
	    	double tempo = sData.get(m).getMaxLogProbability(0);
	    	modelSum += tempo;
	    	modelScores[m] =  tempo;
	    	if(type.equalsIgnoreCase("IL")){
	    		tempo = rsData.get(m).getMaxLogProbability(0);
	    		rmodelSum += tempo;
	    		rmodelScores[m] = tempo;
	    	}
	    }

	    modelScore = (modelSum/(sData.size()-1));
	    if(type.equalsIgnoreCase("IL")){
	    	rmodelScore = (rmodelSum/(sData.size()-1));
	    }
	    if(type.equalsIgnoreCase("IL")){
	    	if(rmodelScore > modelScore)          // choose between forward and reversed for each model
	    	{
	    		modelScore = rmodelScore;
	    		modelScores = rmodelScores;
	    		reversed = 1;
	    	}
	    	else {
	    		reversed = 0;
	    	}
	    }else{
	    	reversed = 0;
	    }

	    double[] groupScores = new double[sData.size()];
	    for(int i=1; i <sData.size(); i++){
	    	groupScores[i-1] = modelScores[i];
	    }

	    double[] deficits = getDeficits(group.Best_Score, groupScores);

	    int[] InteriorMinDist = getMinEditDistance(group,sData,type,reversed,true);
		int[] FullMinDist = getMinEditDistance(group,sData,type,reversed,false);

		boolean[] cutoffs = getCutoffs(group,groupScores,InteriorMinDist);
		double[] cutoffscores = getCutoffScores(group,groupScores,InteriorMinDist);

	    List<SequenceResult> seqRes = new ArrayList<SequenceResult>();
		for(int m = 0; m < sData.size() - 1; m++) {
			SequenceResult seqR;
			if(reversed == 1){
				seqR = new BasicSequenceResult(sData.get(m + 1), groupScores[m], deficits[m], InteriorMinDist[m], FullMinDist[m],reversed,cutoffs[m],cutoffscores[m],rsData.get(m + 1).correspondences);
			}else{
				seqR = new BasicSequenceResult(sData.get(m + 1), groupScores[m], deficits[m], InteriorMinDist[m], FullMinDist[m],reversed,cutoffs[m],cutoffscores[m],sData.get(m + 1).correspondences);
			}
			seqRes.add(seqR);
		}

		LoopResult result = new BasicLoopResult(model_id, reversed, "sig", seqRes);

		for(SequenceResult res : seqRes) {
			res.setLoopResult(result);
		}

	    //Try to clean up the memory we used
	    System.gc();
	    System.runFinalization();
	    System.gc();

	    return seqRes;

	}

	public static List<Sequence> doParse(List<Sequence> sData, String nodeInfo, int range, boolean fullModelText, boolean calculateCorrespondences)
	{
		return doParse(sData,nodeInfo,range,fullModelText,calculateCorrespondences,false);
	}

	// Overloaded doParse that can take model/node data as a string instead of a file name
	// It loops through all Sequence objects in sData.
	// It adds node data just once, to a new Sequence object called S.  Then it adds the organism and letters from each Sequence object.
	// Then it runs the parseSequence method of the Sequence class.
	// It runs through all nodes to store max log probability information
	// Finally, it calls the showParse method of the first node of the model, which accumulates a
	// very, very wide version of the alignments in row-column format.
	// Many columns of the alignment then need to get removed later.

	public static List<Sequence> doParse(List<Sequence> sData, String nodeInfo, int range, boolean fullModelText, boolean calculateCorrespondences, boolean saveNodeProbs) 	{
		Node current;
		List<Double> mProbs = new Vector<Double>();
		double mlp;

		// System.out.println("Alignment.doParse: nodeInfo " + nodeInfo);

		Sequence S = new Sequence("","");                    // blank sequence to use repeatedly
		if (fullModelText){
			S.addNodeDataModelText(nodeInfo);				 // add model data from string
		} else {
			S.addNodeData(nodeInfo);	                     // add model data from file
		}

		sData.get(0).parseData = ((Node)S.first).header(); // add a header line

		Sequence firstS = sData.get(1);          // first sequence, which matches the model
		firstS.setNucleotides();
		firstS.setArrays();

		for (int i = 1; i < sData.size(); i++) {
			S.organism = sData.get(i).organism;             // focus on one sequence
			S.letters  = sData.get(i).letters;

			// don't know why there is an extra * at the end of the sequences
			// System.out.println(String.format("Alignment.doParse: organism %s",S.organism));
			// System.out.println(String.format("Alignment.doParse: letters  %s",S.letters));

			// print the stack to know what called this
			// StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
			// StringBuilder sb = new StringBuilder();
			// for (StackTraceElement element : stackTrace) {
			// 	System.out.println(element.toString());
			// }

			S.setNucleotides();                                    // strip dashes from sequence
			S.setArrays();                                         // define cti, itc, convert letters to numbers

			// System.out.println(String.format("Alignment.doParse S.nucleotides %s", S.nucleotides));

			// keep copies of cti and itc for the first sequence, which must match the model
			S.ctiFirst = new int[firstS.cti.length];
			S.itcFirst = new int[firstS.itc.length];

			for (int j = 0; j < S.ctiFirst.length; j++)
				S.ctiFirst[j] = firstS.cti[j];
			for (int j = 0; j < S.itcFirst.length; j++)
				S.itcFirst[j] = firstS.itc[j];

			mlp = S.parseSequence(range);                             	      // parse this sequence
			sData.get(i).appendProbabilities(mlp); 							  // save mlp

			// We also need to store the parse information somewhere!  All we have is a parse sequence.
			// We need to store the max log probabilities too

			mProbs.clear();
			if(saveNodeProbs)
			{
				current = S.first;
				while(current != null)
				{
					mProbs.add(new Double(current.optimalMaxLogProb));
					//				System.out.println("Alignment.doParse optimalMaxLogProb for a node is "+current.optimalMaxLogProb);
					current =  current.next;
				}
				sData.get(i).appendNodeProbabilities(mProbs);
					System.out.println("doParse: Alignment.doParse actual MLP is "+mlp);
					System.out.println("doParse: mProbs[] is                     "+mProbs);
					System.out.println("doParse: retrieve mProbs gives           "+sData.get(i).getMaxNodeLogProbabilityOf(0,0));
				sData.get(i).parseData = ((InitialNode)S.first).showParse(S.nucleotides);
			}

			if (calculateCorrespondences) {
				String correspondences = ((InitialNode)S.first).showCorrespondences(S.nucleotides);
				correspondences = correspondences.replace("JAR3D_aligns_to", "aligns_to_JAR3D");
				String SF = "Sequence_"+i;
				correspondences = correspondences.replace("SSS",SF);
				sData.get(i).correspondences = correspondences;
				sData.get(i).correspondences += "Sequence_"+i+" has_name "+sData.get(i).organism.replace(" ","_")+"\n";
				sData.get(i).correspondences += "Sequence_"+i+" has_score "+sData.get(i).getMaxLogProbability(0)+"\n";
			}
		}

		if (calculateCorrespondences) {
			for (int i = 1; i < sData.size(); i++)
				{
//					sData.get(i).correspondences += "Sequence_"+i+" has_name "+sData.get(i).organism.replace(" ","_")+"\n";
				}
			for (int i = 1; i < sData.size(); i++)
				{
//					sData.get(i).correspondences += "Sequence_"+i+" has_score "+sData.get(i).getMaxLogProbability(0)+"\n";
				}
		}

		return sData;
	}

	//Overloaded doParse that can take model/node data as a string instead of a file name
	public static List<Sequence> calculateTotalProbability(List<Sequence> sData, String nodeInfo, int range, boolean fullModelText)
	{
        double totalProb;
        int userange;

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

		for (int i = 1; i < sData.size(); i++)              // loop through sequences
		{
			S.organism = sData.get(i).organism;             // focus on one sequence
			S.letters  = sData.get(i).letters;
			S.setNucleotides();                             // strip dashes from sequence
			S.setArrays();                                  // define cti, itc, convert letters to numbers

			// keep copies of cti and itc for the first sequence, which must match the model
			S.ctiFirst = new int[firstS.cti.length];
			S.itcFirst = new int[firstS.itc.length];

			for (int j = 0; j < S.ctiFirst.length; j++)
				S.ctiFirst[j] = firstS.cti[j];
			for (int j = 0; j < S.itcFirst.length; j++)
				S.itcFirst[j] = firstS.itc[j];

			if (range == 0)
				userange = S.letters.length();
			else
				userange = range;

			totalProb = S.calculateTotalProbability(userange); // calculate total probability for this sequence
			sData.get(i).totalProbability = totalProb;      // record total probability for each sequence
		}
		return sData;
	}


	public static int[] getMinEditDistance(MotifGroup group, List<Sequence> sData, String type, int rotation, boolean interior){
		// now works with J3, J4, etc.

		// load sequences from the current motif group
		Vector<Sequence> modsData = Alignment.parseFastaText(group.Sequences,0,0);

		// Calculate edit distances
		// we assume that the sequences in sData are already rotated appropriately
		int[][] EditDistances;
		EditDistances = SimpleAlign.calcEditDistances(sData,modsData,false,interior);

		// for each sequence in sData, return the minimum edit distance to model sequences
		int[] MinDist = new int[EditDistances.length];
		for (int i = 0; i < EditDistances.length; i++) {
			MinDist[i] = ArrayMath.min(EditDistances[i]);
		}
		return MinDist;
	}

	public static boolean[] getCutoffs(MotifGroup group, double[] groupScores, int[] InteriorMinDist){
		boolean[] cutoffs = new boolean[InteriorMinDist.length];
		for(int i = 0; i < InteriorMinDist.length; i++){
			if(groupScores[i]>=group.Cutoffs[0] && InteriorMinDist[i]<=group.Cutoffs[1]
					&& group.Cutoffs[2]*InteriorMinDist[i]-group.Cutoffs[3]*groupScores[i]<=group.Cutoffs[4]){
				cutoffs[i] = Boolean.TRUE;
			}else if(InteriorMinDist[i]==0){
				cutoffs[i] = Boolean.TRUE;
			}else{
				cutoffs[i] = Boolean.FALSE;
			}
		}
		return cutoffs;
	}

	public static double[] getCutoffScores(MotifGroup group, double[] groupScores, int[] InteriorMinDist){
		double[] cutoffscores = new double[InteriorMinDist.length];
		for(int i = 0; i < InteriorMinDist.length; i++){
			cutoffscores[i] = 100 * (group.Cutoffs[2] * InteriorMinDist[i] - group.Cutoffs[3] * groupScores[i] - group.Cutoffs[4]) / -group.Cutoffs[5];
			if(cutoffscores[i] > 100){
				cutoffscores[i] = 100;
			}
			if(cutoffscores[i] > 0 && (groupScores[i] < group.Cutoffs[0] || InteriorMinDist[i] > group.Cutoffs[1])){
				cutoffscores[i] = 0;
			}
		}
		return cutoffscores;
	}

	public static double[] getDeficits(double best_score, double[] groupScores) {
		double [] deficits =  new double[groupScores.length];
		for(int i = 0; i < groupScores.length; i++) {
			deficits[i] = Math.max(0, best_score - groupScores[i]);
		}
		return deficits;
	}
}