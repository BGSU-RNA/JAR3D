import java.io.*;
import java.util.*; 

/**
 * This is the Alignment class, it has methods that are used to facilitate aligning sequences
 * @author meg pirrung
 *
 */
public class Alignment {

	public static Vector JAR3D(String UserDir, String FastaFile, String ModelFile, int numSequences, int DNA, int range) 
	{
		System.setProperty("user.dir",UserDir);
//		 System.out.println(System.getProperty("user.dir"));
		Vector sequenceData = Alignment.loadFastaColumnsDNA(FastaFile,0,0,DNA); 
		sequenceData = Alignment.doParse(sequenceData,numSequences,ModelFile,range);
		Alignment.displayAlignmentFASTA(sequenceData,numSequences);	
		return sequenceData;
	}

	public static Vector loadFasta(String fileName)
	{
		Vector sData = new Vector();
		sData = loadFastaColumns(fileName,0,0);
		return sData;
	}

	/**
	 * This method
	 * @param fileName
	 * @param StartCol
	 * @param EndCol
	 * @return
	 */

	public static Vector loadFastaColumns(String fileName, int StartCol, int EndCol)
	{
		Vector sData = new Vector();
		sData = loadFastaColumnsDNA(fileName,StartCol,EndCol,0);
		return sData;
	}

	
	public static Vector loadFastaColumnsDNA(String fileName, int StartCol, int EndCol, int DNA)
		{
		Vector sData = new Vector();
		String temp="";
		String organism="";
		String letters="";
		sData.add(new Sequence("","")); // save space for header information to be added later as element 0 of this vector
		BufferedReader rdr;
		try {
			String curDir = System.getProperty("user.dir");
			File f1 = new File (curDir + "\\sequences");
			File[] seqfiles = f1.listFiles();
			int i = 0;
			while(!seqfiles[i].getName().equals(fileName))
				i++;

			rdr = new BufferedReader(new FileReader(seqfiles[i]));

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

//			System.out.println(organism);
//			System.out.println(letters);
			
			System.out.println("Read fasta file "+fileName);
			rdr.close();
		}
		catch (IOException e) {
			System.out.println("Could not open fasta file");
			System.out.println(e);
		}
		return sData;
	}

	/**
	 * This method parses using a sequence supplied by a text field
	 * @param fastaText
	 * @param StartCol
	 * @param EndCol
	 * @return
	 */
	public static Vector parseFastaText(String fastaText, int StartCol, int EndCol)
	{
		Vector sData = new Vector();
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
//			System.out.println(organism);
//			System.out.println(letters);
			
			//System.out.println("Read fasta file "+fileName);

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
	 * This method reverses the sequences in sequenceData and returns
	 * @param numSequences
	 * @param sequenceData
	 * @return
	 */
	public static Vector reverse(int numSequences, Vector sequenceData)
	{
		Vector reverseSData = new Vector();

		String temp="";
		String organism="";
		String letters="";
		String second = "";
		String first  = "";
		
		reverseSData.add(new Sequence("","")); // save space for header information to be added later as element 0 of this vector
		
		for (int i=1; i < Math.min(sequenceData.size(),numSequences+1); i++)
		{
			temp = ((Sequence)sequenceData.elementAt(i)).letters;
			organism = ((Sequence)sequenceData.elementAt(i)).organism;
			
			StringTokenizer st = new StringTokenizer(temp,"*");
			second  = st.nextToken();
			first   = st.nextToken();
			letters = first + "*" + second;

			reverseSData.add(new Sequence(organism,letters));
			
//			System.out.println(letters);
		}

		/*
		try {
			String curDir = System.getProperty("user.dir");
			File f1 = new File (curDir + "\\sequences");
			File[] seqfiles = f1.listFiles();
			int i = 0;
			while(!seqfiles[i].getName().equals(fileName))
				i++;

			rdr = new BufferedReader(new FileReader(seqfiles[i]));

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

			sData.add(new Sequence(organism.substring(1,organism.length()),letters));
//			System.out.println(organism);
//			System.out.println(letters);
			
			System.out.println("Read fasta file "+fileName);
			rdr.close();
		}
		catch (IOException e) {
			System.out.println("Could not open fasta file");
			System.out.println(e);
		}
*/
		
		return reverseSData;
	}

	/**
	 * This method initiates parsing
	 * @param sData
	 * @param numSequences
	 * @param nodeFileName
	 * @param range
	 * @return
	 */
	public static Vector doParse(Vector sData, int numSequences, String nodeFileName, int range)
	{
		long start, stop, elapsed;
		Node current;
		Vector mProbs = new Vector();

		Sequence S = new Sequence("","");                        // blank sequence to use repeatedly 
		S.addNodeData(nodeFileName);	                         // only add node data once

		((Sequence)sData.elementAt(0)).parseData = ((InitialNode)S.first).header();

		Sequence firstS = (Sequence)sData.elementAt(1);          // first sequence, which matches the model
		firstS.setNucleotides();
		firstS.setArrays();
		
		for (int i = 1; i < Math.min(sData.size(),numSequences+1); i++) 
		{
			S.organism = ((Sequence)sData.elementAt(i)).organism;             // focus on one sequence
			S.letters  = ((Sequence)sData.elementAt(i)).letters;
			
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

			start = System.currentTimeMillis();
			S.parseSequence(range);                             	      // parse this sequence
			stop = System.currentTimeMillis();
			elapsed = stop - start;

//			System.out.println("Alignment.doParse parsing took "+elapsed/1000+" seconds");
			// We also need to store the parse information somewhere!  All we have is a parse sequence.
			// We need to store the max log probabilities too
			//	pData.add(((InitialNode)S.first).showParse(S.nucleotides) + " " + ((Vector)(((Sequence)sData.elementAt(i)).maxLogProbs.get(0))).get(0));

			mProbs = new Vector();
			current = S.first;
//System.out.println(((InitialNode)S.first).showParse(S.nucleotides)+" ");
			while(current != null)
			{
				mProbs.add(new Double(current.optimalMaxLogProb));
//System.out.println(current.mytype+" "+current.optimalMaxLogProb);
				current =  current.next;
			}
			((Sequence)sData.elementAt(i)).maxLogProbs.add(mProbs);
			if (S.first.optimalMaxLogProb < -99999999)
			{
//				System.out.println("Alignment.doParse: No good parse found, score "+S.first.optimalMaxLogProb);
				current = S.first.next;
				while (current != null)
				{
					if (current.previous.currentMaxLogProb < -99999999 && current.currentMaxLogProb > -999999999)
					{
//						System.out.println("Alignment.doParse: best max Log Prob "+current.previous.currentMaxLogProb+" at "+current.previous.mytype+" "+(current.previous.leftIndex+1)+" "+(current.previous.rightIndex+1));
					}
					current = current.next;				
				}
			}
			
			((Sequence)sData.elementAt(i)).parseData = ((InitialNode)S.first).showParse(S.nucleotides);
			
		}
		return sData;
	}

	/**
	 * This method generates synthetic sequences based on the model selected
	 * @param numSequences is the number of sequences to synthesize
	 * @param modelNum is the number corresponding to which model to use (will probably be a text file soon)
	 * @return
	 */
	public static Vector generateSyntheticData(int numSequences,int range)
	{
		Sequence S = new Sequence("Synthetic organism","AAAA");
		S.addNodeData("");

		Vector sData = new Vector();

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
//			System.out.println(seq);
			sData.add(new Sequence("Synthetic organism",seq));
		}

		return sData;
	}


	/**
	 * This method strips out repetitive dashes used to align sequences
	 * @param pData this is a vector of plain string sequences
	 * @return
	 */
	public static int[] stripDash(Vector pData)
	{
		boolean found = false;
		int[] mask = new int[((String)pData.get(0)).length()];
		for(int i = 0; i<mask.length; i++)
			mask[i] = 0;
		for(int j = 0; j < mask.length; j++)                       // look at each column
		{
			found = false;                                         // has a non-dash been found?

			//System.out.println("Alignment.stripDash Column = "+j);
			//System.out.println("Alignment.stripDash Number of sequences = "+pData.size());

			
			for(int i = 0; i < pData.size() && found != true; i++)
			{
				//System.out.println(" "+i + " " + (String)pData.get(i));
				if(((String)pData.get(i)).charAt(j) != '-')
					found = true;                                  // a non-dash has been found
			}
			//System.out.println();
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
				if(((String)pData.get(i)).charAt(j) != '(')
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
				if(((String)pData.get(i)).charAt(j) != ')')
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
				if(((String)pData.get(i)).charAt(j) != '[')
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
				if(((String)pData.get(i)).charAt(j) != ']')
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
	public static void displayAlignment(Vector sData, int numSequences)
	{
		
		Vector pData = new Vector();
		
//		System.out.println("Alignment.displayAlignment "+sData.size());

		for (int i = 0; i < Math.min(numSequences+1,sData.size()); i++)
		{
			pData.add(((Sequence)sData.get(i)).parseData);
//System.out.println(pData.elementAt(i));
		}
		System.out.println("Displaying alignment ----------------"+pData.size());


		int[] mask = stripDash(pData);

		//System.out.println("Alignment Mask:");
		//for(int i = 0; i < mask.length; i++)
		//	System.out.print(mask[i]);
		//System.out.println();

		System.out.println("Alignment from Java parser:");
		for(int j = 0; j < pData.size(); j++)
		{
			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					System.out.print(((String)pData.get(j)).charAt(i));
			}		
			if(j == 0)
			{
				System.out.println(" "+((Sequence)sData.elementAt(j)).organism);
			}
			else
			{
			System.out.print(" "+((Sequence)sData.elementAt(j)).organism + " ");
			for(int x = 0; x < ((Sequence)sData.elementAt(j)).maxLogProbs.size(); x++)
				System.out.print(((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(x)).get(0)+" ");
			System.out.println();
			}
		}
	}

	/**
	 * This method displays the alignment in FASTA format
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static void displayAlignmentFASTA(Vector sData, int numSequences)
	{
		
		Vector pData = new Vector();
		
		for (int i = 0; i < Math.min(numSequences+1,sData.size()); i++)
		{
			pData.add(((Sequence)sData.get(i)).parseData);
//System.out.println(pData.elementAt(i));
		}
		System.out.println("Displaying alignment ----------------"+pData.size());

		int[] mask = stripDash(pData);

		//System.out.println("Alignment Mask:");
		//for(int i = 0; i < mask.length; i++)
		//	System.out.print(mask[i]);
		//System.out.println();

		System.out.println("Alignment from Java parser:");
		for(int j = 0; j < pData.size(); j++)
		{
			if(j == 0)
			{
				System.out.println("Mask");
			}
			else
			{
				System.out.print(">"+((Sequence)sData.elementAt(j)).organism + " ");
				for(int x = 0; x < ((Sequence)sData.elementAt(j)).maxLogProbs.size(); x++)
					System.out.print(((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(x)).get(0));
				System.out.println();
			}
			for(int i = 0; i < mask.length; i++)
			{
				if (mask[i] == 0)
				{
					char a;
					a = ((String)pData.get(j)).charAt(i);
//					if ((a!='{') && (a!='}') && (a!='<') && (a!='>')) 
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
	public static void displayMaxLogProbs(Vector sData, int numSequences, int Motif, int R)
	{
		
/*		Vector pData = new Vector();
		
		for (int i = 0; i < Math.min(numSequences+1,sData.size()); i++)
		{
			pData.add(((Sequence)sData.get(i)).parseData);
//System.out.println(pData.elementAt(i));
		}
		*/
		
		System.out.println("Displaying alignment ---------------- "+sData.size());

		System.out.println("Alignment from Java parser:");
		for(int j = 1; j < sData.size(); j++)
		{
			System.out.print(Motif+" "+R+" ");
			for(int x = 0; x < ((Sequence)sData.elementAt(j)).maxLogProbs.size(); x++)
				System.out.print(((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(x)).get(0)+" ");
			System.out.print(((Sequence)sData.elementAt(j)).organism + " ");
			System.out.println();
		}

	}

	
	/**
	 * This method displays the alignment
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static Vector getAlignment(Vector sData, int numSequences)
	{
		Vector pData = new Vector();
		
		for (int i = 0; i < Math.min(sData.size(),numSequences+1); i++)
		{
			pData.add(((Sequence)sData.get(i)).parseData);
//			System.out.println((String)pData.get(i));
		}
		
		int[] mask = stripDash(pData);
		String alnm = "";
		Vector alignmentVect = new Vector();

		alnm = "";
		for(int i = 0; i < mask.length; i++)
		{
			if(mask[i] == 0)
				alnm += ((String)pData.get(0)).charAt(i);
		}		
		alignmentVect.add(alnm);                          // paste in header line
		
		for(int j = 1; j < pData.size(); j++)
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
			alignmentVect.add(alnm);
			alignmentVect.add(">"+((Sequence)sData.elementAt(j)).organism);
			alignmentVect.add("Score = "+((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(0)).get(0));
		}
		return alignmentVect;	
	}	
	
	/**
	 * This method displays the alignment
	 * @param pData contains plain string sequences
	 * @param sData contains sequence objects
	 */
	public static Vector getSortedAlignment(Vector seqNames, Vector modNames, int numSequences, int range)
	{
		Vector alignmentVect = new Vector();
		Vector sData = new Vector();
		Vector pData = new Vector();
		double[] modelSums = new double[modNames.size()];
		double[] modelScores = new double[modNames.size()];
		double sum = 0;
		
		// put all sequences in the same Vector
		for(int l = 0; l < seqNames.size(); l++)
		{
			sData.addAll(Alignment.loadFasta((String)seqNames.get(l)));
		}
		
		// parse all sequences against models
		for(int k = 0; k < modNames.size(); k++)
		{
			sData = Alignment.doParse(sData,numSequences,(String)modNames.get(k),30);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < ((Sequence)sData.elementAt(m)).maxLogProbs.size(); x++)
			{
				sum = modelSums[x]; // get current sum for this model (x)
				
				sum = sum + Double.parseDouble((String)((Vector)((Sequence)sData.elementAt(m)).maxLogProbs.get(x)).get(0)); // get score for this sequence(m)
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
			pData.add(((Sequence)sData.get(i)).parseData);

		int[] mask = stripDash(pData);
		String alnm = "";
		
		for(int j = 0; j < pData.size(); j++)
		{
			alnm = "";

			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += ((String)pData.get(j)).charAt(i);
			}		
			alnm += " "+((Sequence)sData.elementAt(j)).organism+" ";
			
			for(int x = 0; x < ((Sequence)sData.elementAt(j)).maxLogProbs.size(); x++)
				alnm += modNames.get(indices[x]) + " score: " + ((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(indices[x])).get(0) + " ";
			alignmentVect.add(alnm);
		}
		return alignmentVect;
	}
	
	public static String getSortedHLAlignment(Vector sData, Vector modNames, int numSequences, int range)
	{
		Vector alignmentVect = new Vector();
		Vector pData = new Vector();
		double[] modelSums = new double[modNames.size()];
		double[] modelScores = new double[modNames.size()];
		Vector shortModNames = new Vector();
		Vector tinyModNames = new Vector();
        String scores = "";

		if(((String)modNames.get(0)).contains("http"))
		{
		// remove http:// from model names
			for(int k = 0; k < modNames.size(); k++)
			{
				shortModNames.add(((String)modNames.get(k)).substring(33,((String)modNames.get(k)).length()-4));
			}
		}
		else
		{
			shortModNames = new Vector(modNames);
		}

		for (int k=0; k< modNames.size(); k++)
		{
			tinyModNames.add(((String)modNames.get(k)).substring(0,8));
		}
		
		// parse all sequences against models
		for(int k = 0; k < modNames.size(); k++)
		{
			System.out.println("Alignment.getSortedHLAlignment: Model is "+modNames.get(k));
			sData  = Alignment.doParse(sData,numSequences,(String)modNames.get(k),range);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < ((Sequence)sData.elementAt(m)).maxLogProbs.size(); x++)
			{
				String temp = String.valueOf(((Vector)((Sequence)sData.elementAt(m)).maxLogProbs.get(x)).get(0)); // get score for this sequence(m)
				double tempo = Double.parseDouble(temp);
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
			scores += fmt.format(" %12.6f", Math.max(modelScores[g],-9999));	
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
			pData.add(((Sequence)sData.get(i)).parseData);

		int[] mask = stripDash(pData);
		String alnm = "";
		
		for(int j = 0; j < pData.size(); j++)
		{
			alnm = "";

			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += ((String)pData.get(j)).charAt(i);
			}		
			alnm += " "+((Sequence)sData.elementAt(j)).organism+" ";
			
			for(int x = 0; x < ((Sequence)sData.elementAt(j)).maxLogProbs.size(); x++)
			{
				fmt = new Formatter();
				fmt.format("%10.6f", ((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(indices[x])).get(0));
				alnm += tinyModNames.get(indices[x]) + " score: " + fmt + " ";
			}
			alignmentVect.add(alnm);
		}
		return scores;
	}	

	public static double[] getSortedILAlignment(Vector sData, Vector modNames, int numSequences, int range)
	{
		Vector alignmentVect = new Vector();                   // alignment lines to output
		Vector pData = new Vector();                           // parse data
		double[] modelSums = new double[modNames.size()];      // sum of alignment scores
		double[] rmodelSums = new double[modNames.size()];     // sum with sequences reversed
		double[] modelScores = new double[modNames.size()];
		double[] rmodelScores = new double[modNames.size()];
		Vector shortModNames = new Vector();                   // for easier display
		Vector tinyModNames = new Vector();                    // for even easier display
		int[] reversed = new int[modNames.size()];             // is best model reversed?
 //       String scores = "";
        double[] scores = new double[2*modNames.size()];       // all scores computed

        Vector rsData = Alignment.reverse(numSequences, sData);  // reversed sequence data
		
		if(((String)modNames.get(0)).contains("http"))         // look online for models
		{
		// remove http:// from model names
			for(int k = 0; k < modNames.size(); k++)
			{
				shortModNames.add(((String)modNames.get(k)).substring(33,((String)modNames.get(k)).length()-4));
			}
		}
		else
		{
			shortModNames = new Vector(modNames);
		}

		for (int k=0; k< modNames.size(); k++)
		{
			tinyModNames.add(((String)modNames.get(k)).substring(0,8));
		}
		
		// parse sequence data in sData against models
		for(int k = 0; k < modNames.size(); k++)
		{
			// System.out.println("Alignment.getSortedILAlignment " + modNames.get(k));
			sData  = Alignment.doParse(sData,numSequences,(String)modNames.get(k),range);
			// System.out.println("Alignment.getSortedILAlignment " + modNames.get(k));
			rsData = Alignment.doParse(rsData,numSequences,(String)modNames.get(k),range);
		}

		// add up model scores for each sequence
		for(int m = 0; m < sData.size(); m++)
		{
			for(int x = 0; x < ((Sequence)sData.elementAt(m)).maxLogProbs.size(); x++)
			{
				// there is no good reason for having to do these crazy manipulations
				// in order to get the value of a double variable, but at least this works
				String temp = String.valueOf(((Vector)((Sequence)sData.elementAt(m)).maxLogProbs.get(x)).get(0)); // get score for this sequence(m)
				double tempo = Double.parseDouble(temp);   
				// System.out.print(tempo+"  ");
				modelSums[x] += tempo;

				temp = String.valueOf(((Vector)((Sequence)rsData.elementAt(m)).maxLogProbs.get(x)).get(0)); // get score for this sequence(m)
				tempo = Double.parseDouble(temp);
				rmodelSums[x] += tempo;
			}
			// System.out.println("");
		}
		
		Formatter fmt = new Formatter(); 
		
		System.out.println("Average score for each model unsorted: ");
		for(int g = 0; g < modelSums.length; g++)
		{
			modelScores[g] = (modelSums[g]/(sData.size()-1));
			rmodelScores[g] = (rmodelSums[g]/(sData.size()-1));
			System.out.format("%s          %12.6f     ",tinyModNames.get(g),modelScores[g]);
			System.out.format("%s reversed %12.6f\n",tinyModNames.get(g),rmodelScores[g]);
			
			fmt = new Formatter();
			scores[2*g]   = Math.max(modelScores[g],-9999);
			scores[2*g+1] = Math.max(rmodelScores[g],-9999);

//			fmt = new Formatter();
//			scores += fmt.format(" %12.6f", Math.max(rmodelScores[g],-9999));
			
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
        	sData  = Alignment.doParse(sData,numSequences,(String)modNames.get(indices[0]),range);
        	System.out.println("parsing forward again "+(String)modNames.get(indices[0])+" "+indices[0]);
    		for (int i = 0; i < sData.size(); i++)
    			pData.add(((Sequence)sData.get(i)).parseData);  
        }
        else {
        	rsData = Alignment.doParse(rsData,numSequences,(String)modNames.get(indices[0]),range);
        	System.out.println("parsing reversed again "+(String)modNames.get(indices[0])+" "+indices[0]);
    		for (int i = 0; i < sData.size(); i++)
    			pData.add(((Sequence)rsData.get(i)).parseData); 
        }
		
		

		int[] mask = stripDash(pData);
		String alnm = "";

		for(int j = 0; j < Math.min(pData.size(),30); j++)    // at most 30 sequences from each
		{
			alnm = "";

			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += ((String)pData.get(j)).charAt(i);
			}		
			
			if (j > 0) {
			alnm += " "+((Sequence)sData.elementAt(j)).organism+" ";
			for(int x = 0; x < modNames.size(); x++)
			{
				fmt = new Formatter();
				if (reversed[indices[x]] == 0)
					fmt.format("%12.6f", ((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(indices[x])).get(0));
				else
					fmt.format("%12.6f", ((Vector)((Sequence)rsData.elementAt(j)).maxLogProbs.get(indices[x])).get(0));
				alnm += tinyModNames.get(indices[x]) + " score: " + fmt + " ";
			}
			}
        	System.out.println(alnm);
		}

		return scores;
	}	
	
	public static void makeHTMLAlignment(Vector sData, int numSequences)
	{
		Vector pData = new Vector();
		String html = "";
		
		html+= "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"\n";
	    html+= "\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n";
	    html+= "<html xmlns=\"http://www.w3.org/1999/xhtml\">\n";
		html+= "<head>\n";
		html+= "<H1><font face=\"Tahoma\">Alignment</font></H1>\n";
		html+= "</head>\n";
		html+= "<body bgcolor=\"#ccffff\">\n";
		html+= "<title>Alignment</title>\n";
		html+= "<table bgcolor=\"#ffffff\" border=\"0\" cellpadding=3>\n";
		html+= "<tr bgcolor=\"#0066cc\">\n";
		html+= "<th><font face=\"Tahoma\" color=\"#ffffff\">Alignment</font></th>\n";
		html+= "<th><font face=\"Tahoma\" color=\"#ffffff\">Organism Info</font></th>\n";
		html+= "<th><font face=\"Tahoma\" color=\"#ffffff\">Score</font></th>\n";
		html+= "</tr>\n";
		
		for (int i = 0; i < Math.min(sData.size(),numSequences+1); i++)
			pData.add(((Sequence)sData.get(i)).parseData);
	
		int[] mask = stripDash(pData);
		String alnm = "";
		Vector alignmentVect = new Vector();
	
		for(int j = 0; j < Math.min(pData.size(),numSequences+1); j++)
		{
			alnm = "";
			for(int i = 0; i < mask.length; i++)
			{
				if(mask[i] == 0)
					alnm += ((String)pData.get(j)).charAt(i);
			}		
			if(j%2 == 0)
			{
				html+= "<tr bgcolor=\"#ebebeb\">\n";
			}
			else
				html+= "<tr bgcolor=\"#ffffff\">\n";
			
			html+= "<td><font size=2 face=\"Courier\">\n";
			html+= alnm;
			html+= "</font></td>\n";
			html+= "<td><font size=2 face=\"Tahoma\">\n";
			html+= ((Sequence)sData.elementAt(j)).organism;
			html+= "</font></td>\n";
			html+= "<td><font size=2 face=\"Tahoma\">\n";
			
			for(int x = 0; x < ((Sequence)sData.elementAt(j)).maxLogProbs.size(); x++)
				html+= ((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(x)).get(0);
			
			html+= "</font></td>\n";
			html+="</tr>\n";
			
		}
		html+="</table>\n";
		try {
	        BufferedWriter out = new BufferedWriter(new FileWriter("alignment.html"));
	        out.write(html);
	        out.close();
	    } catch (IOException e) {
	    }
	}

	public static void printAlignment(Vector aData, int numChars)
	{
		for(int f = 0; f < ((String)aData.get(1)).length()/numChars+1; f++)
		{
			for(int j = 0; j < aData.size(); j++)
			{
				if((f+1)*numChars > ((String)aData.get(j)).length())
				{
					System.out.println(((String)aData.get(j)).substring(f*numChars, ((String)aData.get(j)).length()));
				}
				else
					System.out.println(((String)aData.get(j)).substring(f*numChars, (f+1)*numChars));
			}
			System.out.println();
		}
	}

	public static String tempdisplayAlignmentFASTA(Vector sData, int numSequences)
	{
		String temp = "";
		Vector pData = new Vector();
		
		for (int i = 0; i < Math.min(numSequences+1,sData.size()); i++)
		{
			pData.add(((Sequence)sData.get(i)).parseData);
//System.out.println(pData.elementAt(i));
		}
		temp += "Displaying alignment ----------------"+pData.size()+"\n";


		int[] mask = stripDash(pData);

		//System.out.println("Alignment Mask:");
		//for(int i = 0; i < mask.length; i++)
		//	System.out.print(mask[i]);
		//System.out.println();

		System.out.println("Alignment from Java parser:");
		for(int j = 0; j < pData.size(); j++)
		{
			if(j == 0)
			{
				temp += "Mask\n";
			}
			else
			{
				temp += ">"+((Sequence)sData.elementAt(j)).organism + " ";
				for(int x = 0; x < ((Sequence)sData.elementAt(j)).maxLogProbs.size(); x++)
					temp += ((Vector)((Sequence)sData.elementAt(j)).maxLogProbs.get(x)).get(0);
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

}