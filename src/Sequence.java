import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.net.*;
import java.io.*;

/**
 * Sequence class of the Parser project, Used to store the branching data structure
 * of different Node types
 * @author meg pirrung
 *
 */
public class Sequence {
	long start, stop, elapsed;
	String organism;
	String letters;            // as read from the fasta file
	String nucleotides = "";   // with gaps stripped out
	int ctiFirst[];            // maps fasta column to sequence index for reference sequence
	int itcFirst[];
	int cti[];                 // maps fasta column to sequence index
	int itc[];                 // maps sequence index to fasta column

	Node first, last;          // first and last nodes of parsing model

	int code[];                // code for letters A, C, G, U
	Vector maxLogProbs;
    String parseData;
    String correspondences;
    
	/**
	 * Constructor for Sequence class
	 * sets up a new sequence using parameters read from the fasta file
	 * @param org holds the organism name
	 * @param let holds the sequence of characters as read from the file
	 */
	public Sequence(String org, String let) {
		organism = org;
		letters = let;
		maxLogProbs = new Vector();

		/*                                        This will be done just before parsing
		setNucleotides();
		setArrays();

		ctiFirst = new int[cti.length];
		itcFirst = new int[itc.length];

		for (int i = 0; i < cti.length; i++)
			ctiFirst[i] = cti[i];
		for (int i = 0; i < itc.length; i++)
			itcFirst[i] = itc[i];
		 */
	}

	/**
	 * Constructor for Sequence class that adds cti and itc information from first sequence
	 * sets up a new sequence using parameters read from the fasta file
	 * @param org holds the organism name
	 * @param let holds the sequence of characters as read from the file
	 */
	public Sequence(String org, String let, Sequence firstSeq) {
		organism = org;
		letters = let;
		setNucleotides();
		setArrays();
		maxLogProbs = new Vector();

		ctiFirst = new int[firstSeq.cti.length];
		itcFirst = new int[firstSeq.itc.length];

		for (int i = 0; i < firstSeq.cti.length; i++)
			ctiFirst[i] = firstSeq.cti[i];
		for (int i = 0; i < firstSeq.itc.length; i++)
			itcFirst[i] = firstSeq.itc[i];

	}

	/**
	 * Copy constructor
	 * @param aSeq
	 */
	public Sequence(Sequence aSeq)
	{
		this(aSeq.organism, aSeq.letters);
		setNucleotides();
		setArrays();
	}

	/**
	 * This method is used to strip dashes from the fasta file for
	 * a plain nucleotide string
	 */
	void setNucleotides()
	{
		nucleotides = "";
		StringTokenizer st = new StringTokenizer(letters,"-");
		while(st.hasMoreTokens())
		{
			nucleotides+=st.nextToken();
		}	 
	}

	/**
	 * This method sets up the different arrays used in the Sequence class
	 * cti is the column (of fasta file) to index (in 3D structure) array
	 * itc is the index to column array
	 * code is an array of ints used to hold numerical codes for each nucleotide letter
	 */
	void setArrays()
	{
		cti = new int[letters.length()];
		int[] tem = new int[letters.length()];
		int i = 0;
		for(int c = 0; c < letters.length(); c++)
		{
			cti[c] = i;                     // gaps and non-gaps get mapped to current index
			if(letters.charAt(c)!= '-')
			{
				tem[i] = c;                 // map i back to the column
				i++;                        // increase i of non-gap letters
			}
		}

		itc = new int[i];                   // make itc have the right length

		for (int k = 0; k < i; k++)
			itc[k] = tem[k];                // fill in the values found earlier

		/*
		System.out.println("Sequence.setArrays");

		for(int c = 0; c < letters.length(); c++)
		{
			System.out.print(letters.charAt(c));
		}

		System.out.println();

		System.out.print("cti: ");
		for(int c = 0; c < letters.length(); c++)
		{
			System.out.print(cti[c]+" ");
		}

		System.out.println();

		System.out.print("itc: ");
		for(int c = 0; c < itc.length; c++)
		{
			System.out.print(itc[c]+" ");
		}

		System.out.println();
		 */

		code = new int[nucleotides.length()];

		for(int c = 0; c < nucleotides.length(); c++)
		{
			switch(nucleotides.charAt(c))
			{
			case 'A':
				code[c] = 0;
				break;
			case 'C':
				code[c] = 1;
				break;
			case 'G':
				code[c] = 2;
				break;
			case 'U':
				code[c] = 3;
				break;
			case 'a':
				code[c] = 0;
				break;
			case 'c':
				code[c] = 1;
				break;
			case 'g':
				code[c] = 2;
				break;
			case 'u':
				code[c] = 3;
				break;
			case 'T':
				code[c] = 3;
				break;
			case 't':
				code[c] = 3;
				break;
			case '*':
				code[c] = 4;
				break;

				// other characters get the default code 0, which means they are converted to A's
			}
		}
	}

	/**
	 * This method is used to set up the sequence of nodes in a linked list type structure that
	 * also has a branching structure, depending on whether the list is traversed using next
	 * or child nodes.
	 * It takes the text of a model as input
	 * @param modelFileName a string to hold which model file is to be read
	 */
	void addNodeDataModelText(String modelText)
	{
		start = System.currentTimeMillis();
		Node current = new Node();
		Stack branchingNodes = new Stack(); // keep track of current junction
		int newchild = 0;
		int charCnt = 0;
		String nodeLine = "";
		String dataString = "";
		String comment = "";
		int numFix = 0;
		StringTokenizer element;
		boolean flag = true;
		boolean isHairpin = true;
		double[] nData;
		long start1 = System.currentTimeMillis();
		int nodeNumber = 1;
		int lineNum = 0;
				
		String[] modelArray;
		modelArray = modelText.split("\\n");
					
		while(lineNum < modelArray.length)
		{
			nodeLine = modelArray[lineNum];
			nodeLine = nodeLine.replace(" ","");
			
			if (!nodeLine.equals("//") && (!nodeLine.equals("")))
			{
				//System.out.println(nodeLine);
				//Vector numData; // to hold data straight from txt file
				Vector numDatas = new Vector(); // to hold arrays of doubles
				Vector stringData = new Vector();
				if(nodeLine.indexOf("//") != -1)
				{
					comment = nodeLine.substring(nodeLine.indexOf("//"),nodeLine.length());
					nodeLine = nodeLine.substring(0,nodeLine.indexOf("//"));
				}
				StringTokenizer st = new StringTokenizer(nodeLine,"|");
				dataString = st.nextToken(); // the name of the node
				char type  = dataString.charAt(0);
				char iType = dataString.charAt(2);
				//System.out.print("Read node string: " + nodeLine);
				//System.out.println("Node type: " + type);
				//System.out.println("Comments: " + comment);

				// allow comment lines in the data file, and blank lines too

				if(iType == 'a' ) // character definition
				{
					dataString = st.nextToken();
					element = new StringTokenizer(dataString, ",");
					charCnt = element.countTokens();
//						System.out.println("Number of characters to parse with: " + charCnt);
					char[] codesFromDef = new char[charCnt];
					for(int i = 0; i < charCnt; i++)
						codesFromDef[i] = element.nextToken().charAt(0);
				}
				else
				{
					while(st.hasMoreTokens())
					{
						dataString = st.nextToken(); // each individual set of data
						//System.out.println(dataString);
						dataString = dataString.replace("[","*");
						dataString = dataString.replace("]","*");
						dataString = dataString.replace(",","*");
						dataString = dataString.replace(" ","");

						element = new StringTokenizer(dataString, "*");

						nData = new double[element.countTokens()-1];
						element.nextToken();
						for(int i = 0; i < nData.length; i++)
						{
							String nT;
							nT = element.nextToken();
							//System.out.println(nT);
							nData[i] = Double.parseDouble(nT);
						}
						numDatas.add(nData);	
					}
					
					switch(type)
					{	
					case 'I':
						switch(iType)
						{
						case 'i':
							current = new InitialNode(current, (double[])numDatas.get(0), (double[])numDatas.get(1), (double[])numDatas.get(2), (double[])numDatas.get(3),(int)((double[])numDatas.get(4))[0],(int)((double[])numDatas.get(5))[0]);
							current.number = nodeNumber;
							nodeNumber++;
							if (flag)
							{
								first = current;
								flag = false;
							}
							break;
						case 't':
							double[][] interactionParams = new double[charCnt][charCnt];
							for(int i = 0; i < charCnt; i++)
							{
								for(int j = 0; j < charCnt; j++)
								{
									interactionParams[i][j] = ((double[])numDatas.get(1))[charCnt*i + j];
								}
							}
							if(isHairpin)
							{
								((HairpinNode)current).addInteraction((int)((double[])numDatas.get(0))[0],(int)((double[])numDatas.get(0))[1], interactionParams);
							}
							else
								((ClusterNode)current).addInteraction((int)((double[])numDatas.get(0))[0],(int)((double[])numDatas.get(0))[1], interactionParams);
							break;
						case 's':
							if(isHairpin)
							{
								((HairpinNode)current).addInsertion((int)((double[])numDatas.get(0))[0],((double[])numDatas.get(1)),((double[])numDatas.get(2)));
							}
							else
								((ClusterNode)current).addInsertion((int)((double[])numDatas.get(0))[0],((double[])numDatas.get(1)),((double[])numDatas.get(2)));
							break;
						default:
						}
						break;
					case 'B':
						double[][] pairProb = new double[charCnt][charCnt];
						for(int i = 0; i < charCnt; i++)
						{
							for(int j = 0; j < charCnt; j++)
							{
								pairProb[i][j] = ((double[])numDatas.get(1))[charCnt*i + j];
							}
						}
						current = new BasepairNode(current, ((double[])numDatas.get(0))[0], pairProb, (double[])numDatas.get(2), (double[])numDatas.get(3), (double[])numDatas.get(4), (double[])numDatas.get(5), (int)((double[])numDatas.get(6))[0],(int)((double[])numDatas.get(7))[0]);
						current.number = nodeNumber;
						nodeNumber++;
						break;
					case 'F':
						current = new FixedNode(current, ((double[])numDatas.get(0))[0], (double[])numDatas.get(1), (int)((double[])numDatas.get(2))[0], (int)((double[])numDatas.get(3))[0],(int)((double[])numDatas.get(4))[0]);
						current.number = nodeNumber;
						nodeNumber++;
						break;
					case 'J':
						current = new JunctionNode(current, (int)((double[])numDatas.get(0))[0], (int)((double[])numDatas.get(1))[0],(int)((double[])numDatas.get(2))[0]);
						current.number = nodeNumber;
						nodeNumber++;
						break;
					case 'H':
						current = new HairpinNode(current, (int)((double[])numDatas.get(0))[0], (int)((double[])numDatas.get(1))[0], (int)((double[])numDatas.get(2))[0], ((double[])numDatas.get(3))[0]);
						current.number = nodeNumber;
						nodeNumber++;
						isHairpin = true;
						break;
					case 'C':
						current = new ClusterNode(current, ((double[])numDatas.get(0))[0], (int)((double[])numDatas.get(1))[0], (int)((double[])numDatas.get(2))[0], (int)((double[])numDatas.get(3))[0], (int)((double[])numDatas.get(4))[0], ((double[])numDatas.get(5))[0]);
						current.number = nodeNumber;
						nodeNumber++;
						isHairpin = false;
						break;
					case 'A':
						int numAlternatives = (int)((double[])numDatas.get(0))[0];
						double[] priorProb = new double[numAlternatives];
						for(int i = 0; i < numAlternatives; i++)
						{
							priorProb[i] = ((double[])numDatas.get(1))[i];
						}
						current = new AlternativeNode(current,numAlternatives,priorProb,(int)((double[])numDatas.get(2))[0],(int)((double[])numDatas.get(3))[0]);
						current.number = nodeNumber;
						nodeNumber++;
					default:
					}
				}
				lineNum++;

			}
			else
				lineNum++;
		}

		stop = System.currentTimeMillis();
		elapsed = stop - start1;
//		System.out.println("Reading model data file time: " + elapsed + " milliseconds");

		start1 = System.currentTimeMillis();


		// Common to all models we might make

		last = current;
		last.next = null; //The last element does not have a next, 
		// so it is set to null so that it can be
		// be used as a stopping case


		current = last; // The current is set to last because you
		// work from the back to assign previous nodes

		while(current.previous != null) // go until the end of the list
		{
			current.previous.next = current; // assign the previous node's next
			current = current.previous;		 // current equals the one before
		}

		stop = System.currentTimeMillis();
		elapsed = stop - start1;
//		System.out.println("Connecting nodes 1 time: " + elapsed + " milliseconds");
		start1 = System.currentTimeMillis();


		current = first; // set current equal to first for traversal and
		// assignment of children

		/*
		while(current != null) // Go until end of list
		{
			if(current.getType().equals("ClusterNode"))
				((ClusterNode)current).normalize();
			if(current.getType().equals("HairpinNode"))
				((HairpinNode)current).normalize();
			current = current.next;
		}
		 */

		while(current != null) // Go until end of list
		{
			if(current.getType().equals("ClusterNode"))
				((ClusterNode)current).useFileNormalize();
//				((ClusterNode)current).normalize();
			if(current.getType().equals("HairpinNode"))
				((HairpinNode)current).useFileNormalize();
//				((HairpinNode)current).normalize();
			current = current.next;
		}


		stop = System.currentTimeMillis();
		elapsed = stop - start1;
//		System.out.println("Connecting nodes 2 time: " + elapsed + " milliseconds");
		start1 = System.currentTimeMillis();



		current = first; // set current equal to first for traversal and
		// assignment of children

		while(current != null) // Go until end of list
		{
			if(current.getType().equals("JunctionNode") || current.getType().equals("AlternativeNode")) // if the type
			{	
				branchingNodes.push(current); // add this junction to the stack
				((BranchingNode)current).children.add(current.next); // add the next
			}					             						// as a child
			else if(current.getType().equals("HairpinNode")||current.getType().equals("eNode"))
			{
				if(current.getType().equals("eNode"))
					((AlternativeNode)branchingNodes.peek()).eNodes.add(current);

				newchild = 0;
				while(!(branchingNodes.isEmpty()) && newchild == 0)
				{
					((BranchingNode)branchingNodes.peek()).branches--;

					if(((BranchingNode)branchingNodes.peek()).branches == 0)
					{
						if(current.getType().equals("AlternativeNode"))
						{
							for(int i = 0; i <((AlternativeNode)branchingNodes.peek()).eNodes.size();i++)
								((BranchingNode)branchingNodes.peek()).children.add(((eNode)((AlternativeNode)branchingNodes.peek()).eNodes.get(i)).next);
						}
						branchingNodes.pop();
						//System.out.println("Done with one brachingNode.");
					}
					else
					{
						// add the one after the current(a hairpin or eNode)
						((BranchingNode)branchingNodes.peek()).children.add(current.next);
						newchild = 1;
					}
				}
			}
			else 
			{
				current.child = current.next;
				// System.out.println(current.getType());
			}// end else if
			current = current.next;
		}//end



		// For debugging purposes, show what junction goes with what children
		/*
		current = first; 

		while(current != null) // Go until end of list
		{
			int lI, rI;
			if(current.getType().equals("JunctionNode"))
			{
//				if ((current.leftIndex != ((Node)((BranchingNode)current).children.get(0)).leftIndex) || (current.rightIndex != ((Node)((BranchingNode)current).children.get(1)).rightIndex))
				{
					System.out.println(current.leftIndex+" "+current.rightIndex);
					for(int i = 0; i < ((BranchingNode)current).children.size(); i++)
					{
						System.out.print("Child left index: " + ((Node)((BranchingNode)current).children.get(i)).leftIndex);
						System.out.println(" Child right index: " + ((Node)((BranchingNode)current).children.get(i)).rightIndex);
					}

				}
			}
			current = current.next;
		}
		 */

		stop = System.currentTimeMillis();
		elapsed = stop - start1;
//		System.out.println("Connecting nodes 3 time: " + elapsed + " milliseconds");

		elapsed = stop - start;
//		System.out.println("Node data added time: " + elapsed + " milliseconds");
	}// end method addNodedataModelText


	/**
	 * This method reads a model file into a text string and then calls addNodeDataModelText
	 * @param modelFileName a string to hold which model file is to be read
	 */
	void addNodeData(String modelFileName)
	{
		String modelText = "";
		String nodeLine = "";

		BufferedReader rdr;
		
		try {
			if(!(modelFileName.contains("http")))
			{
				String curDir = System.getProperty("user.dir");
		        curDir = curDir.replace(File.separator + "bin","");

				try
				{
					File f2 = new File(curDir + File.separator + "Models" + File.separator + modelFileName);
					rdr = new BufferedReader(new FileReader(f2));
				}
				catch(FileNotFoundException e)
				{
					System.out.println("Reading model file from absolute path");
					File f3 = new File(modelFileName);
					rdr = new BufferedReader(new FileReader(f3));
//					System.out.println("Sequence.addNodeData: Model file not found.");
//					System.exit(0);
				}
			}
			else
			{
				// if the model filename has http in it, use URL stuff, if not, look for a local file.
				System.out.println("Working online...");
				URL dirurl = new URL("http://rna.bgsu.edu/JAR3D/");
		        URLConnection dircon = dirurl.openConnection();
		        rdr = new BufferedReader(new InputStreamReader(dircon.getInputStream()));
	
		        String ln = rdr.readLine();
				while(ln != null)
					{
					ln = rdr.readLine();
					}
	
				dirurl = new URL(modelFileName);
		        dircon = dirurl.openConnection();
		        rdr = new BufferedReader(new InputStreamReader(dircon.getInputStream())); 
			}
			
			nodeLine = rdr.readLine();
			nodeLine = nodeLine.replace(" ","");
			
			while(nodeLine != null)
			{
				nodeLine = nodeLine.replace(" ","");
				modelText = modelText + nodeLine + "\n";
				nodeLine = rdr.readLine();
			}		
			
			addNodeDataModelText(modelText);
			}
		catch (IOException e) {
			System.out.println("Sequence.addNodeData: Could not open model file");
			System.out.println(e);
		}
	}
	
	
	/**
	 * This method is used to set up the sequence of nodes in a linked list type structure that
	 * also has a branching structure, depending on whether the list is traversed using next
	 * or child nodes
	 * @param modelFileName a string to hold which model file is to be read
	 */
	void addNodeDataOld(String modelFileName)
	{
		start = System.currentTimeMillis();
		Node current = new Node();
		Stack branchingNodes = new Stack(); // keep track of current junction
		int newchild = 0;
		int charCnt = 0;
		String nodeLine = "";
		String dataString = "";
		String comment = "";
		int numFix = 0;
		StringTokenizer element;
		boolean flag = true;
		boolean isHairpin = true;
		double[] nData;
		long start1 = System.currentTimeMillis();

		BufferedReader rdr;
		
		try {
			if(!(modelFileName.contains("http")))
			{
				String curDir = System.getProperty("user.dir");
		        curDir = curDir.replace(File.separator + "bin","");
//				System.out.println("Sequence.addNodeData: Current directory is "+curDir);

//				File f1 = new File (curDir + "\\..\\models");
//				File[] modfiles = f1.listFiles();
				int z = 0;
	
				try
				{
//					while(!modfiles[z].getName().equals(modelFileName))
//						z++;
				}
				catch(ArrayIndexOutOfBoundsException e)
				{
					System.out.println("Sequence.addNodeData: Model file not found.");
					System.exit(0);
				}
				File f2 = new File(curDir + File.separator + "models" + File.separator + modelFileName);
				rdr = new BufferedReader(new FileReader(f2));
			}
			else
			{
				// if the model filename has http in it, use URL stuff, if not, look for a local file.
				System.out.println("Working online...");
				URL dirurl = new URL("http://rna.bgsu.edu/JAR3D/");
		        URLConnection dircon = dirurl.openConnection();
		        rdr = new BufferedReader(new InputStreamReader(dircon.getInputStream()));
	
		        String ln = rdr.readLine();
				while(ln != null)
					{
					ln = rdr.readLine();
					}
	
				dirurl = new URL(modelFileName);
		        dircon = dirurl.openConnection();
		        rdr = new BufferedReader(new InputStreamReader(dircon.getInputStream())); 
			}
			
			nodeLine = rdr.readLine();
			nodeLine = nodeLine.replace(" ","");
			
			while(nodeLine != null)
			{
				if (!nodeLine.equals("//") && (!nodeLine.equals("")))
				{
					//System.out.println(nodeLine);
					//Vector numData; // to hold data straight from txt file
					Vector numDatas = new Vector(); // to hold arrays of doubles
					Vector stringData = new Vector();
					if(nodeLine.indexOf("//") != -1)
					{
						comment = nodeLine.substring(nodeLine.indexOf("//"),nodeLine.length());
						nodeLine = nodeLine.substring(0,nodeLine.indexOf("//"));
					}
					StringTokenizer st = new StringTokenizer(nodeLine,"|");
					dataString = st.nextToken(); // the name of the node
					char type  = dataString.charAt(0);
					char iType = dataString.charAt(2);
					//System.out.print("Read node string: " + nodeLine);
					//System.out.println("Node type: " + type);
					//System.out.println("Comments: " + comment);

					// allow comment lines in the data file, and blank lines too

					if(iType == 'a' ) // character definition
					{
						dataString = st.nextToken();
						element = new StringTokenizer(dataString, ",");
						charCnt = element.countTokens();
//						System.out.println("Number of characters to parse with: " + charCnt);
						char[] codesFromDef = new char[charCnt];
						for(int i = 0; i < charCnt; i++)
							codesFromDef[i] = element.nextToken().charAt(0);
					}
					else
					{
						while(st.hasMoreTokens())
						{
							dataString = st.nextToken(); // each individual set of data
							//System.out.println(dataString);
							dataString = dataString.replace("[","*");
							dataString = dataString.replace("]","*");
							dataString = dataString.replace(",","*");
							dataString = dataString.replace(" ","");

							element = new StringTokenizer(dataString, "*");

							nData = new double[element.countTokens()-1];
							element.nextToken();
							for(int i = 0; i < nData.length; i++)
							{
								String nT;
								nT = element.nextToken();
								//System.out.println(nT);
								nData[i] = Double.parseDouble(nT);
							}
							numDatas.add(nData);	
						}
						
						switch(type)
						{	
						case 'I':
							switch(iType)
							{
							case 'i':
								current = new InitialNode(current, (double[])numDatas.get(0), (double[])numDatas.get(1), (double[])numDatas.get(2), (double[])numDatas.get(3),(int)((double[])numDatas.get(4))[0],(int)((double[])numDatas.get(5))[0]);
								if (flag)
								{
									first = current;
									flag = false;
								}
								break;
							case 't':
								double[][] interactionParams = new double[charCnt][charCnt];
								for(int i = 0; i < charCnt; i++)
								{
									for(int j = 0; j < charCnt; j++)
									{
										interactionParams[i][j] = ((double[])numDatas.get(1))[charCnt*i + j];
									}
								}
								if(isHairpin)
								{
									((HairpinNode)current).addInteraction((int)((double[])numDatas.get(0))[0],(int)((double[])numDatas.get(0))[1], interactionParams);
								}
								else
									((ClusterNode)current).addInteraction((int)((double[])numDatas.get(0))[0],(int)((double[])numDatas.get(0))[1], interactionParams);
								break;
							case 's':
								if(isHairpin)
								{
									((HairpinNode)current).addInsertion((int)((double[])numDatas.get(0))[0],((double[])numDatas.get(1)),((double[])numDatas.get(2)));
								}
								else
									((ClusterNode)current).addInsertion((int)((double[])numDatas.get(0))[0],((double[])numDatas.get(1)),((double[])numDatas.get(2)));
								break;
							default:
							}
							break;
						case 'B':
							double[][] pairProb = new double[charCnt][charCnt];
							for(int i = 0; i < charCnt; i++)
							{
								for(int j = 0; j < charCnt; j++)
								{
									pairProb[i][j] = ((double[])numDatas.get(1))[charCnt*i + j];
								}
							}
							current = new BasepairNode(current, ((double[])numDatas.get(0))[0], pairProb, (double[])numDatas.get(2), (double[])numDatas.get(3), (double[])numDatas.get(4), (double[])numDatas.get(5), (int)((double[])numDatas.get(6))[0],(int)((double[])numDatas.get(7))[0]);
							break;
						case 'J':
							current = new JunctionNode(current, (int)((double[])numDatas.get(0))[0], (int)((double[])numDatas.get(1))[0],(int)((double[])numDatas.get(2))[0]);
							break;
						case 'H':
							current = new HairpinNode(current, (int)((double[])numDatas.get(0))[0], (int)((double[])numDatas.get(1))[0], (int)((double[])numDatas.get(2))[0], ((double[])numDatas.get(3))[0]);
							isHairpin = true;
							break;
						case 'C':
							current = new ClusterNode(current, ((double[])numDatas.get(0))[0], (int)((double[])numDatas.get(1))[0], (int)((double[])numDatas.get(2))[0], (int)((double[])numDatas.get(3))[0], (int)((double[])numDatas.get(4))[0], ((double[])numDatas.get(5))[0]);
							isHairpin = false;
							break;
						case 'A':
							int numAlternatives = (int)((double[])numDatas.get(0))[0];
							double[] priorProb = new double[numAlternatives];
							for(int i = 0; i < numAlternatives; i++)
							{
								priorProb[i] = ((double[])numDatas.get(1))[i];
							}
							current = new AlternativeNode(current,numAlternatives,priorProb,(int)((double[])numDatas.get(2))[0],(int)((double[])numDatas.get(3))[0]);
						default:
						}
					}
					nodeLine = rdr.readLine();

				}
				else
					nodeLine = rdr.readLine();
			}
			//rdr.close();
			rdr.close();
		}
		catch (IOException e) {
			System.out.println("Sequence.addNodeData: Could not open model file");
			System.out.println(e);
		}

		stop = System.currentTimeMillis();
		elapsed = stop - start1;
//		System.out.println("Reading model data file time: " + elapsed + " milliseconds");

		start1 = System.currentTimeMillis();


		// Common to all models we might make

		last = current;
		last.next = null; //The last element does not have a next, 
		// so it is set to null so that it can be
		// be used as a stopping case


		current = last; // The current is set to last because you
		// work from the back to assign previous nodes

		while(current.previous != null) // go until the end of the list
		{
			current.previous.next = current; // assign the previous node's next
			current = current.previous;		 // current equals the one before
		}

		stop = System.currentTimeMillis();
		elapsed = stop - start1;
//		System.out.println("Connecting nodes 1 time: " + elapsed + " milliseconds");
		start1 = System.currentTimeMillis();


		current = first; // set current equal to first for traversal and
		// assignment of children

		/*
		while(current != null) // Go until end of list
		{
			if(current.getType().equals("ClusterNode"))
				((ClusterNode)current).normalize();
			if(current.getType().equals("HairpinNode"))
				((HairpinNode)current).normalize();
			current = current.next;
		}
		 */

		while(current != null) // Go until end of list
		{
			if(current.getType().equals("ClusterNode"))
				((ClusterNode)current).useFileNormalize();
//				((ClusterNode)current).normalize();
			if(current.getType().equals("HairpinNode"))
				((HairpinNode)current).useFileNormalize();
//				((HairpinNode)current).normalize();
			current = current.next;
		}


		stop = System.currentTimeMillis();
		elapsed = stop - start1;
//		System.out.println("Connecting nodes 2 time: " + elapsed + " milliseconds");
		start1 = System.currentTimeMillis();



		current = first; // set current equal to first for traversal and
		// assignment of children

		while(current != null) // Go until end of list
		{
			if(current.getType().equals("JunctionNode") || current.getType().equals("AlternativeNode")) // if the type
			{	
				branchingNodes.push(current); // add this junction to the stack
				((BranchingNode)current).children.add(current.next); // add the next
			}					             						// as a child
			else if(current.getType().equals("HairpinNode")||current.getType().equals("eNode"))
			{
				if(current.getType().equals("eNode"))
					((AlternativeNode)branchingNodes.peek()).eNodes.add(current);

				newchild = 0;
				while(!(branchingNodes.isEmpty()) && newchild == 0)
				{
					((BranchingNode)branchingNodes.peek()).branches--;

					if(((BranchingNode)branchingNodes.peek()).branches == 0)
					{
						if(current.getType().equals("AlternativeNode"))
						{
							for(int i = 0; i <((AlternativeNode)branchingNodes.peek()).eNodes.size();i++)
								((BranchingNode)branchingNodes.peek()).children.add(((eNode)((AlternativeNode)branchingNodes.peek()).eNodes.get(i)).next);
						}
						branchingNodes.pop();
						//System.out.println("Done with one brachingNode.");
					}
					else
					{
						// add the one after the current(a hairpin or eNode)
						((BranchingNode)branchingNodes.peek()).children.add(current.next);
						newchild = 1;
					}
				}
			}
			else 
			{
				current.child = current.next;
				// System.out.println(current.getType());
			}// end else if
			current = current.next;
		}//end



		// For debugging purposes, show what junction goes with what children
		/*
		current = first; 

		while(current != null) // Go until end of list
		{
			int lI, rI;
			if(current.getType().equals("JunctionNode"))
			{
//				if ((current.leftIndex != ((Node)((BranchingNode)current).children.get(0)).leftIndex) || (current.rightIndex != ((Node)((BranchingNode)current).children.get(1)).rightIndex))
				{
					System.out.println(current.leftIndex+" "+current.rightIndex);
					for(int i = 0; i < ((BranchingNode)current).children.size(); i++)
					{
						System.out.print("Child left index: " + ((Node)((BranchingNode)current).children.get(i)).leftIndex);
						System.out.println(" Child right index: " + ((Node)((BranchingNode)current).children.get(i)).rightIndex);
					}

				}
			}
			current = current.next;
		}
		 */






		stop = System.currentTimeMillis();
		elapsed = stop - start1;
//		System.out.println("Connecting nodes 3 time: " + elapsed + " milliseconds");

		elapsed = stop - start;
//		System.out.println("Node data added time: " + elapsed + " milliseconds");
	}// end method addNodedata

	/**
	 * Prints out nodes and their data
	 */
	void printNodeData()
	{
		// print out information about all nodes
		Node current;

		current = first; // set current to first for traversal

		System.out.println("Printing first to last");

		while(current != null)
		{
			System.out.println("Type: " + current.getType() + " Params: "+ current.getParams());
			if(current.getType().equals("JunctionNode") || current.getType().equals("AlternativeNode"))
			{
				for(int i = 0; i < ((BranchingNode)current).children.size(); i++)
				{
					System.out.println("Child params: " + ((Node)((BranchingNode)current).children.get(i)).getParams());
				}
			}			
			current = current.next;
		}


		System.out.println();
		System.out.println("Printing last to first");

		// print node data, last node to first

		/*		current = last; // set current to last for backwards traversal

		while(current != null)
		{
			System.out.println("Type: " + current.getType() + " Params: "+ current.getParams());
			if(current.getType().equals("JunctionNode"))
			{
				for(int i = 0; i < ((JunctionNode)current).children.size(); i++)
				{
					System.out.println("Child params: " + ((Node)((JunctionNode)current).children.get(i)).getParams());
				}
			}					
			current = current.previous;
		}
		 */
	}

	/**
	 * Calls generate of the first node in the list and returns a string of the resulting
	 * generated characters
	 * @param del this is a boolean used by certain nodes to determine deletion probabilities
	 * @return
	 */
	String generate(boolean del)
	{
		return first.generate(true); // true is a placeholder
		// how to decide where to start off
		// variable length helices
	}


	/**
	 * This method converts the sequence into an array of numbers [0,1,2,3]
	 * and allocates space in each node for maxprob, according to the length 
	 * of this sequence
	 */
	void parseSequence(int range)
	{
		start = System.currentTimeMillis();

		Node current = last; 
		while(current != null) // go until the end of the list
		{
			current.currentMaxLogProb = -1d/0d;
			current.setMinMax(range, this);
			current = current.previous;		 // current equals the one before
		}

//		System.out.println("Sequence.parseSequence: nucleotides 0 to "+(nucleotides.length()-1)+" "+nucleotides);
		
		// CYK algorithm for determining maximum log probability parse
		for (int p=1; p <= nucleotides.length(); p++)         // length of subsequence
		{
			long start2 = System.currentTimeMillis();
			for(int i=0; i+p <= nucleotides.length(); i++)   // starting point of subsequence
			{
				int j = i+p-1;                          // right side of this subsequence
				current = last;					        // start with last node
				while (current != null)
				{
					current.computeMaxLogProb(this,i,j);         // determine how this node would generate subsequence i to j
					current = current.previous;
				}// while
			}// for		     	
			stop = System.currentTimeMillis();
			elapsed = stop - start2;
//			System.out.println("Subsequences of length "+p+" took "+elapsed+" milliseconds");

		}// for	
		stop = System.currentTimeMillis();
		elapsed = stop - start;
//		System.out.println("Time to compute max log prob: " + elapsed+ " milliseconds");


		start = System.currentTimeMillis();
		
		((InitialNode)first).traceback(0,nucleotides.length()-1);
		stop = System.currentTimeMillis();
		elapsed = stop - start;
//		System.out.println("Time to traceback: " + elapsed+ " milliseconds");
	}// end parseSequence()
	
	/**
	 * This method looks in the local folder for HL, IL, or JL model names,
	 * but if it can't look there, it looks online at rna.bgsu.edu/JAR3D.
	 */
	
	public static Vector getModelNames(String loopType)
	{
		Vector modelNames = new Vector();
		String listName = loopType + "_Models.txt";

	try {
		String curDir = System.getProperty("user.dir");
        curDir = curDir.replace(File.separator + "bin","");

//		System.out.println("Sequence.getModelNames: Current directory is "+curDir);
//        File f1 = new File (curDir + "\\models");
//	    File[] modfiles = f1.listFiles();

		BufferedReader rdr;

		File f2 = new File(curDir + File.separator + "models" + File.separator + listName);
		rdr = new BufferedReader(new FileReader(f2));

		String fileLine = "";
		fileLine = rdr.readLine();
//		System.out.println(fileLine);
		
		while(fileLine != null)
		{
//			System.out.println(fileLine);
			modelNames.add(fileLine);
			fileLine = rdr.readLine();
		}

		}
	catch (Exception e) {
		try{
			System.out.println("Couldn't find model files locally, looking online.");
			URL dirurl = new URL("http://rna.bgsu.edu/JAR3D/models/"+listName);
	        URLConnection dircon = dirurl.openConnection();
	        BufferedReader rdr = new BufferedReader(new InputStreamReader(dircon.getInputStream()));
	        String ln = rdr.readLine();
			while(ln != null)
				{
					modelNames.add("http://rna.bgsu.edu/JAR3D/models/" + ln);
					ln = rdr.readLine();
				}
		}
		catch (Exception ex)
		{
			System.out.println(ex.getMessage());
		}
	}
	return modelNames;
	}
	
	//Overloaded getModelNames for the new file system
	public static Vector getModelNames(String folder, String modelType, boolean Structured)
	{
		Vector modelNames = new Vector();
		char fsep = File.separatorChar;
		String listName;
		if(Structured){
			listName = folder + fsep + modelType + "_models" + fsep + "structured.txt";
		}else{
			listName = folder + fsep + modelType + "_models" + fsep + "all.txt";
		}

	try {
		String curDir = System.getProperty("user.dir");
        curDir = curDir.replace(File.separator + "bin","");

//		System.out.println("Sequence.getModelNames: Current directory is "+curDir);
//        File f1 = new File (curDir + "\\models");
//	    File[] modfiles = f1.listFiles();

		BufferedReader rdr;

		File f2 = new File(listName);
		rdr = new BufferedReader(new FileReader(f2));

		String fileLine = "";
		fileLine = rdr.readLine();
//		System.out.println(fileLine);
		
		while(fileLine != null)
		{
//			System.out.println(fileLine);
			//remove _model.txt from file name to get group name
			fileLine = fileLine.substring(0,fileLine.length()-10);
			modelNames.add(fileLine);
			fileLine = rdr.readLine();
		}

		}
	catch (Exception e) {
		try{
			System.out.println("Couldn't find model files locally, looking online.");
			URL dirurl = new URL("http://rna.bgsu.edu/JAR3D/models/"+listName);
	        URLConnection dircon = dirurl.openConnection();
	        BufferedReader rdr = new BufferedReader(new InputStreamReader(dircon.getInputStream()));
	        String ln = rdr.readLine();
			while(ln != null)
				{
					modelNames.add("http://rna.bgsu.edu/JAR3D/models/" + ln);
					ln = rdr.readLine();
				}
		}
		catch (Exception ex)
		{
			System.out.println(ex.getMessage());
		}
	}
	return modelNames;
	}


	/**
	 * This method reads a text file and returns the lines as a vector
	 * 
	 */
	
	public static Vector readTextFile(String fileName)
	{
		Vector lineValues = new Vector();

	try {
		String curDir = System.getProperty("user.dir");
        curDir = curDir.replace(File.separator + "bin","");

		BufferedReader rdr;

		File f2 = new File(curDir + File.separator + fileName);
		rdr = new BufferedReader(new FileReader(f2));

		String fileLine = "";
		fileLine = rdr.readLine();
		
		while(fileLine != null)
		{
			lineValues.add(fileLine);
			fileLine = rdr.readLine();
		}

		}
	catch (Exception e) {
			System.out.println(e.getMessage());
		}

	return lineValues;
	}

}// end class

