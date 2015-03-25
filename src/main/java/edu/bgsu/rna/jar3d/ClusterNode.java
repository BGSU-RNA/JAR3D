package edu.bgsu.rna.jar3d;
import java.util.*;

/**
 * This class is used to simulate cluster interactions
 * @author meg pirrung
 *
 */
public class ClusterNode extends BasicNode {

	LinkedList interactions = new LinkedList();
	Vector<InsertionDistribution> insertions;
	int numLeftInt, numRightInt;
	double normalizationZ;
	int[] maxLengths;
	
	//	 this is the probability that this ClusterNode does not generate anything
    double deleteProb;
    double deleteLogProb;
    double notDeleteLogProb;
   
	/**
	 * This is the constructor for ClusterNodes
	 * @param prev contains a pointer to the previous node
	 * @param dProb is the deletion probability
	 * @param numLIn is the number of insertions on the left side
	 * @param numRIn is the number of insertions on the right side
	 */
	public ClusterNode(Node prev, double dProb, int numLIn, int numRIn, int lI, int rI, double nC) {
		super(prev, "ClusterNode", lI, rI);
		numLeftInt = numLIn;
		numRightInt = numRIn;
		normalizationZ = nC;
		deleteProb = dProb;                                   // deletion probability
		deleteLogProb = Math.log(deleteProb);                 // log probability
		notDeleteLogProb = Math.log(1-deleteProb);            // probability of not being deleted
		
		insertions = new Vector<InsertionDistribution>(numLeftInt+numRightInt);
		for (int i = 0; i<numLeftInt+numRightInt; i++)
		{
            // This dummy setting is below changed when an insertion is added
			insertions.add(new InsertionDistribution(new double[]{1.0}, new double[]{0.0, 0.3,0.2,0.1,0.5}));
		}
	}
	
	/**
	 * This method adds an interaction between two bases of the ClusterNode
	 * @param first is the first base
	 * @param second is the second base
	 * @param subsProbs is a matrix of substituion probabilities for the bases
	 */
	void addInteraction(int first,int second,double[][] subsProbs)
	{
		interactions.add(new ClusterInteraction(first-1,second-1,subsProbs));
	}
	
	/**
	 * This method adds insertions between two bases of a ClusterNode
	 * @param ins is the position at which to add the insertion
	 * @param lDist is the length distribution matrix of the insertion
	 * @param letDist is the letter distribution matrix of the insertion
	 */
	void addInsertion(int ins, double[] lDist, double[] letDist)
	{
		if ((ins >= 1) && ((ins < numLeftInt) || (ins >= numLeftInt+1)) && (ins < numLeftInt+numRightInt))
			insertions.set(ins-1, new InsertionDistribution(lDist,letDist));
	}


	/**
	 * This method scores the given sequence based on the insertions and interactions present
	 */
	public void normalize()
	{
		double z = 0;
		double p = 1;
		int numBases = numLeftInt + numRightInt;
		int[] codes = new int[numBases];
		
		System.out.println("ClusterNode.Normalize "+numBases);

		for (int i = 0; i < numBases; i++)
			codes[i] = 0;
				
		maxLengths = new int[insertions.size()]; // maximum possible numbers of insertions
		for(int l = 0; l < insertions.size(); l++)
		{
			maxLengths[l] = ((InsertionDistribution)insertions.get(l)).lengthDist.length-1;
		}

		for(int i = 0; i < (Math.pow(4,numBases)); i++)
		{
			// score codes[] according to the various interactions

			p = 1;
			for(int j = 0; j < interactions.size(); j++)
			{
				p *= ((ClusterInteraction)interactions.get(j)).getSubstProb(codes);
			}
			z += p;     // running total of probabilities
			// push a 1 into codes and carry
			codes[numBases-1] += 1;
			
			int place = numBases-1;
			while((codes[place] > 3) && (place > 0))// carry
			{
				codes[place] = 0;
				place--;
				codes[place]++;
			}
			
		}
		normalizationZ = z;
		System.out.println("ClusterNode.Normalize normalized a node: " + z);
	}

	public void useFileNormalize()
	{
		double z = 0;
		double p = 1;
		int numBases = numLeftInt + numRightInt;
		int[] codes = new int[numBases];
		
		for (int i = 0; i < numBases; i++)
			codes[i] = 0;
				
		maxLengths = new int[insertions.size()]; // maximum possible numbers of insertions
		for(int l = 0; l < insertions.size(); l++)
		{
			maxLengths[l] = ((InsertionDistribution)insertions.get(l)).lengthDist.length-1;
		}
		
//		normalizationZ = 1;
	}

	/**
	 * This method generates interacting bases for the ClusterNode
	 * @return
	 */
	public String[] generateInteractingBases()
	{
		double z = 0;
		double p = 1;

		double limit = Math.random() * normalizationZ;
		
		int numBases = numLeftInt + numRightInt;
		int[] codes = new int[numBases];
		
		for (int i = 0; i < numBases; i++)
			codes[i] = 0;
		
		codes[numBases-1] -= 1;
		
		int i = 0;
		
		while((z <= limit) && (i < Math.pow(4,numBases)))
		{
			// push a 1 into codes and carry
			codes[numBases-1] += 1;

			int place = numBases-1;
			while((codes[place] > 3) && (place > 0))// carry
			{
				codes[place] = 0;
				place--;
				codes[place]++;
			}

			// score codes[] according to the various interactions
			p = 1;
			for(int j = 0; j < interactions.size(); j++)
			{
				p *= ((ClusterInteraction)interactions.get(j)).getSubstProb(codes);
			}
			z += p;
			
			i++;
			
		}
		
		char letters[] = new char[]{'A','C','G','U'};
		
		
// 		take first numLeftInt and return separately from numRightInt
		String leftStr = "";
		for(int f = 0; f < numLeftInt; f++)
			leftStr += letters[codes[f]];
		String rightStr = "";
		for(int g = numLeftInt; g < numRightInt+numLeftInt; g++)
			rightStr += letters[codes[g]];
		
		return (new String[]{leftStr,rightStr});
	}
	
	/**
	 * This method generates insertions for the ClusterNode
	 * @return
	 */
	public Vector<String> generateInsertions()
	{
		Vector<String> ins = new Vector<String>(numLeftInt+numRightInt);
		String temp = "";
		for(int i = 0; i < ins.capacity(); i++)
		{
			 temp = insertions.get(i).generate();
			 ins.add(temp);
		}
		return ins;
	}
	
	public String generate()
	{
		String[] bases = generateInteractingBases();
		
		Vector<String> subseq = generateInsertions();
		
		String left = "";
		String right = "";
		for(int i = 0; i<bases[0].length(); i++)
			left += bases[0].charAt(i)+(String)subseq.get(i);
		for(int i = 0; i<bases[1].length(); i++)
			right += bases[1].charAt(i) + (String)subseq.get(i+numLeftInt);

		System.out.println("ClusterNode: "+left+"----"+right);
		
		return "{" + left + child.generate(true) + right + "}"; // using true as a holder
	}	
	
	
	public void computeMaxLogProb(Sequence seq, int i, int j)
	{
		/**
		 * If the subsequence from i to j is in the range to be parsed by this node, 
		 * Note, however, that i+numLeftInt might exceed j-numRightInt; in that case, this node really must be deleted and generate nothing
		 */
		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax) && (i <= j))
		{
			double p;									// maximum log probability found so far
			double pnew; 								// probability of current insertion possibility
			int a = 0; 									// number of insertions on the left
			int b = 0; 									// number of insertions on the right
			int opta = 0, optb = 0;                     // optimal numbers of insertions on left and right
			int numBases = numLeftInt + numRightInt;    // number of interacting bases in this Cluster
			int[] intCodes = new int[numBases];            	// codes of interacting bases
			int[] insLengths = new int[insertions.size()]; 	// one combination of insertion numbers 
			int[] optInsLengths = new int[insertions.size()]; // optimal combination of insertions
			
			
			// start off with no insertions anywhere, score this possibility
			for(int l = 0; l < insertions.size(); l++)
			{
                insLengths[l] = 0;
			}
			/*
			 * loop through insertion combinations.
			 * pull out interacting base codes
			 * for each one, add up log probabilities for:
			 * 1. the interacting bases, following this model:
			 * 			p = 0;
			 * */
			
			boolean Deleted;
			
			// consider the possibility that this node generates nothing at all
			
			p = deleteLogProb + super.child.getMaxLogProb(i,j); 
			Deleted = true;                          // default, can be changed later

			while(insLengths[0] <= maxLengths[0])
			{
				// if there is enough space between i and j for the int and ins bases,
				if (i+numLeftInt+a <= j-numRightInt-b)                  // used to be < and -b-1, changed 9-29-07 CLZ
				{
					// pull out the numBases interacting base codes
	
					pnew = notDeleteLogProb;
	
					int k = i;                                       // first base of the subsequence being parsed
					for(int m = 0; m < numLeftInt; m++)              // go through insertions on the left
					{
						intCodes[m] = seq.code[k];                   // record code of interacting base
						int[] insCodes = new int[insLengths[m]];     // place for codes of insertions
						for (int q = 0; q < insLengths[m]; q++)      // loop through number of insertions at site m
							insCodes[q] = seq.code[k+1+q];           // record codes of right number of inserted bases here
						pnew += ((InsertionDistribution)insertions.get(m)).computeLogProb(insCodes);
						k += insLengths[m] + 1;
					}
						
					k = j;                                           // last base of the subsequence being parsed
					for(int m = numBases-1; m >= numLeftInt; m--)    // consider insertions on the right
					{
						intCodes[m] = seq.code[k];                   // record codes of interacting bases, starting from the right
						int[] insCodes = new int[insLengths[m-1]];   // (-1)
						for (int q = 0; q < insLengths[m-1]; q++)    // (-1) loop through number of insertions at site m, moving from the *right*
						{
							insCodes[q] = seq.code[k-1-q];           // record code of base here
						}
						pnew += ((InsertionDistribution)insertions.get(m-1)).computeLogProb(insCodes);
						k = k - insLengths[m-1] - 1;                 // (-1)

					}
										
					for(int m = 0; m < interactions.size(); m++)
					{
//						System.out.println("ClusterNode:"+intCodes);
						pnew += ((ClusterInteraction)interactions.get(m)).getLogSubstProb(intCodes);
					}
					
					
					pnew += super.child.getMaxLogProb(i+numLeftInt+a, j-numRightInt-b);
					
					pnew -= Math.log(normalizationZ); // normalize the probability

					// compare to current best possibility and keep if better
					if(pnew > p)
					{
						//should we keep the codes array for insertions?
						p = pnew;
						for (int z = 0; z < insLengths.length; z++)
						{
							optInsLengths[z] = insLengths[z];
						}
						opta = a;
						optb = b;
						Deleted = false;
					}
				}
					
				// push a 1 into insLengths and carry
				// this is how it manages to consider all possible insertion numbers and locations
				// "carrying" is done when the number of insertions at a particular location exceeds the allowed number there
				insLengths[numBases-1]++;

				int place = numBases-1;
				while((insLengths[place] > maxLengths[place]) && (place > 0))// carry
				{
					insLengths[place] = 0;
					place--;
					insLengths[place]++;
				}

				// update the total number of left and right insertions
				a = 0;
				for(int m = 0; m < numLeftInt; m++)
				{
					a += insLengths[m];
				}
					
				b = 0;
				for(int m = numBases-1; m >= numLeftInt; m--)
				{
					b += insLengths[m];
				}

			}
			// for loop that sets maxLogProb[i-super.iMin][j-super.jMin]
			
			maxLogProb[i-iMin][j-jMin] = p;

  			if (Deleted)
  			{
  	  			myGen[i-iMin][j-jMin] = new genData(true);
  			}
  			else
  			{
  				myGen[i-iMin][j-jMin] = new genData(optInsLengths,opta,optb);
  			}
		} // end if loop
	}// end method computMaxLogProb

	
	
	public void traceback(int i, int j)
  	{
  		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j))
  		{
  			setOptimalAndRelease(i,j);
  			
   			if(optimalGen1.deleted)
  			{
  				super.child.traceback(i,j);
  			}
  			else
  			{	  			
  				int numLeftIns = optimalGen1.numLeftIns;
  				int numRightIns = optimalGen1.numRightIns;
  			
  				if ((i+numLeftInt+numLeftIns) > (j-numRightInt-numRightIns))
  				{
  					System.out.println("ClusterNode.traceback:  Child gets "+(i+numLeftInt+numLeftIns)+" to "+(j-numRightInt-numRightIns));
  					System.out.println("i "+i+" j "+j+" numLeftIns "+numLeftIns+" numRightIns "+numRightIns+" numLeftInt "+numLeftInt+" numRightInt "+numRightInt);
  				}

  				super.child.traceback(i+numLeftInt+numLeftIns,j-numRightInt-numRightIns);
  			}
  		}
  		else
  		{
  			System.out.println("ClusterNode.traceback: Cluster node out of range");
  			System.out.println("Indices "+leftIndex+" "+rightIndex);
  		}
  	}

	
	public String showParse(String n)
  	{
		int[] insert = optimalGen1.insLengths;
		int i = optimalGen1.i;
		int j = optimalGen1.j;
				
  		if (optimalGen1.deleted)
  		{
  			String left = "-";
			
			for(int l = 0; l < numLeftInt-1; l++)
				{
					for(int f = 0; f < maxLengths[l]; f++)
					{
						left += "-";
					}
	                left += "-";
				}
			
			String right = "-";
			
			for(int l = numLeftInt+numRightInt-2; l >= numLeftInt; l--)
				{
					for(int f = 0; f < maxLengths[l]; f++)
						right = "-" + right;
	                right = "-" + right;
				}
			return "{" + left + super.child.showParse(n) + right + "}";
  		}
  		else
  		{
			
			String left = n.substring(i,i+1);
	
			for(int l = 0; l < numLeftInt-1; l++)
				{
					left += n.substring(i+1,i+1+insert[l]);
					for(int f = 0; f < maxLengths[l] - insert[l] ; f++)
					{
						left += "-";
					}
					i += insert[l]+1;
	                left += n.substring(i,i+1);
				}
			
			String right = n.substring(j,j+1);
			
			for(int l = numLeftInt+numRightInt-2; l >= numLeftInt; l--)
				{
					right = n.substring(j-(insert[l]),j) + right;
					for(int f = 0; f < maxLengths[l] - insert[l]; f++)
						right = "-" + right;
					j = j - insert[l] - 1;
	                right = n.substring(j,j+1) + right;
				}
			return "{" + left + super.child.showParse(n) + right + "}";
  		}
		
  	}
	
	public String header()
	{		
		String left = "";

		for(int l = 0; l < numLeftInt; l++)
			{
				left += "I";
				for(int f = 0; f < maxLengths[l]; f++)
					left += "-";
			}
		
		String right = "";
		
		for(int l = numLeftInt+numRightInt-2; l >= numLeftInt-1; l--)
			{
				right = "I" + right;
				for(int f = 0; f < maxLengths[l]; f++)
					right = "-" + right;
			}
		return "{" + left + super.child.header() + right + "}";
	}

	public String showCorrespondences(String letters)
  	{
		int[] insert = optimalGen1.insLengths;
		int i = optimalGen1.i;
		int j = optimalGen1.j;
		
  		if (optimalGen1.deleted)
  		{
			return super.child.showCorrespondences(letters);
  		}
  		else
  		{
  			String left = "SSS_Position_" + (i+1) + "_" + letters.charAt(i) + " JAR3D_aligns_to " + "MMM_Node_" + number + "_Position_1" + "\n";

			for(int l = 0; l < numLeftInt-1; l++)
				{
					for (int k = i+1; k <= i+insert[l]; k++)
		  				left += "SSS_Position_" + (k+1) + "_" + letters.charAt(k) + " JAR3D_aligns_to " + "MMM_Node_" + number + "_Position_" + (l+1) + "-" + (l+2) + "_Insertion" + "\n";
					i += insert[l]+1;
					left += "SSS_Position_" + (i+1) + "_" + letters.charAt(i) + " JAR3D_aligns_to " + "MMM_Node_" + number + "_Position_" + (l+2) + "\n";
				}

			String right = "SSS_Position_" + (j+1) + "_" + letters.charAt(j) + " JAR3D_aligns_to " + "MMM_Node_" + number + "_Position_" + (numLeftInt+numRightInt) + "\n";
			
			for(int l = numLeftInt+numRightInt-2; l >= numLeftInt; l--)
				{
					for (int k = j - 1; k >= j - insert[l]; k--)
						right = "SSS_Position_" + (k+1) + "_" + letters.charAt(k) + " JAR3D_aligns_to " + "MMM_Node_" + number + "_Position_" + (l+1) + "-" + (l+2) + "_Insertion" + "\n" + right;
					j = j - insert[l] - 1;
					right = "SSS_Position_" + (j+1) + "_" + letters.charAt(j) + " JAR3D_aligns_to " + "MMM_Node_" + number + "_Position_" + (l+1) + "\n" + right;
				}
 
  			return left + super.child.showCorrespondences(letters) + right;
  		}
  	}


	public void computeTotalProbability(Sequence seq, int i, int j)
	{
		/**
		 * If the subsequence from i to j is in the range to be parsed by this node, 
		 * Note, however, that i+numLeftInt might exceed j-numRightInt; in that case, this node really must be deleted and generate nothing
		 */
		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax) && (i <= j))
		{
			double p;									// total probability found so far
			double pnew; 								// probability of current insertion possibility
			int a = 0; 									// number of insertions on the left
			int b = 0; 									// number of insertions on the right
			int numBases = numLeftInt + numRightInt;    // number of interacting bases in this Cluster
			int[] intCodes = new int[numBases];            	// codes of interacting bases
			int[] insLengths = new int[insertions.size()]; 	// one combination of insertion numbers 
			
			p = 0;      // start with zero probability
			
//			System.out.println("ClusterNode1: " + p);

			// start off with no insertions anywhere, score this possibility
			for(int l = 0; l < insertions.size(); l++)
			{
                insLengths[l] = 0;
			}
			/*
			 * loop through insertion combinations.
			 * pull out interacting base codes
			 * for each one, add up log probabilities for:
			 * 1. the interacting bases, following this model:
			 * 			p = 0;
			 * */

			while(insLengths[0] <= maxLengths[0])
			{
				// if there is enough space between i and j for the interacting and inserted bases,
				if (i+numLeftInt+a <= j-numRightInt-b)               // 
				{
					// pull out the numBases interacting base codes

					
					pnew = 1;

					int k = i;                                       // first base of the subsequence being parsed
					for(int m = 0; m < numLeftInt; m++)              // go through insertions on the left
					{
						intCodes[m] = seq.code[k];                   // record code of interacting base
						int[] insCodes = new int[insLengths[m]];     // place for codes of insertions
						for (int q = 0; q < insLengths[m]; q++)      // loop through number of insertions at site m
							insCodes[q] = seq.code[k+1+q];           // record codes of right number of inserted bases here
						pnew *= ((InsertionDistribution)insertions.get(m)).computeProb(insCodes);
						k += insLengths[m] + 1;
					}

					k = j;                                           // last base of the subsequence being parsed
					for(int m = numBases-1; m >= numLeftInt; m--)    // consider insertions on the right
					{
						intCodes[m] = seq.code[k];                   // record codes of interacting bases, starting from the right
						int[] insCodes = new int[insLengths[m-1]];   // (-1)
						for (int q = 0; q < insLengths[m-1]; q++)    // (-1) loop through number of insertions at site m, moving from the *right*
						{
							insCodes[q] = seq.code[k-1-q];           // record code of base here
						}
						pnew *= ((InsertionDistribution)insertions.get(m-1)).computeProb(insCodes);
						k = k - insLengths[m-1] - 1;                 // (-1)
					}


					for(int m = 0; m < interactions.size(); m++)
					{
						pnew *= ((ClusterInteraction)interactions.get(m)).getSubstProb(intCodes);
					}

					pnew *= super.child.getTotalProb(i+numLeftInt+a, j-numRightInt-b);

					pnew = pnew / normalizationZ;                    // normalize the probability

					p += pnew;
				}
					
				// push a 1 into insLengths and carry
				// this is how it manages to consider all possible insertion numbers and locations
				// "carrying" is done when the number of insertions at a particular location exceeds the allowed number there
				insLengths[numBases-1]++;

				int place = numBases-1;
				while((insLengths[place] > maxLengths[place]) && (place > 0))    // carry
				{
					insLengths[place] = 0;
					place--;
					insLengths[place]++;
				}

				// update the total number of left and right insertions
				a = 0;
				for(int m = 0; m < numLeftInt; m++)
				{
					a += insLengths[m];
				}
					
				b = 0;
				for(int m = numBases-1; m >= numLeftInt; m--)
				{
					b += insLengths[m];
				}

			}
			// for loop that sets maxLogProb[i-super.iMin][j-super.jMin]

			p *= (1-deleteProb);
			
			// consider the possibility that this node generates nothing at all

			p += deleteProb * super.child.getTotalProb(i,j); 

  			totalProb[i-iMin][j-jMin] = p;

		} // end if loop
	}// end method computeTotalProbability


}
