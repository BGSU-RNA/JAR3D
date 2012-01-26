import java.util.LinkedList;
import java.util.Vector;

/**
 * This class is used to simulate Hairpins
 * @author meg pirrung
 *
 */
public class HairpinNode extends BasicNode {

	LinkedList interactions = new LinkedList();
	Vector insertions;
	int numFixed;
	double normalizationZ;
	int[] maxLengths; // maximum possible numbers of insertions
	
	/**
	 * This is the HairpinNode constructor
	 * @param prev contains a pointer to the previous node
	 * @param numfix is the number of fixed bases
	 */
	public HairpinNode(Node prev, int numfix, int lI, int rI) {
		super(prev, "HairpinNode", lI, rI);

		if (numfix > 0)
		{			
			numFixed    = numfix;
//			System.out.println("HairpinNode "+numfix);
			insertions  = new Vector(numFixed);
			for (int i = 0; i<numFixed; i++)
			{
				insertions.add(new InsertionDistribution(new double[]{1.0}, new double[]{0.4, 0.3,0.2,0.1}));
			}
		}
		else
		{
			numFixed = 4;
			insertions  = new Vector(numFixed);
			for (int i = 0; i < numFixed; i++)
			{
				insertions.add(new InsertionDistribution(new double[]{1.0}, new double[]{0.4, 0.3,0.2,0.1}));
				interactions.add(new ClusterInteraction(i,i,new double[][]{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,1}}));
			}
		}
	}
	
	/**
	 * This method adds an interaction between two bases in the hairpin
	 * @param first is the first interacting base
	 * @param second is the second interacting base
	 * @param subsProbs is a matrix of subsitution probabilities for the bases
	 */
	void addInteraction(int first,int second,double[][] subsProbs)
	{
		interactions.add(new ClusterInteraction(first-1,second-1,subsProbs));
	}
	
	/**
	 * This method adds an insertion at the specified position in the hairpin
	 * @param ins is the position of the insertion
	 * @param lDist is the length distribution matrix of the insertion
	 * @param letDist is the letter distribution matrix of the insertion
	 */
	void addInsertion(int ins, double[] lDist, double[] letDist)
	{
		if (ins < numFixed)
			insertions.set(ins-1, new InsertionDistribution(lDist,letDist));
	}

	/**
	 * This method scores the given sequence based on the insertions and interactions present
	 */
	public void normalize()
	{
		double z = 0;
		double p = 1;
		int[] codes = new int[numFixed];
		
		System.out.println("HairpinNode.Normalize "+numFixed);

		maxLengths = new int[insertions.size()]; // maximum possible numbers of insertions
		for(int l = 0; l < insertions.size(); l++)
		{
			maxLengths[l] = ((InsertionDistribution)insertions.get(l)).lengthDist.length-1;
		}

		
		for (int i = 0; i < numFixed; i++)
			codes[i] = 0;
				
		for(int i = 0; i < (Math.pow(((InsertionDistribution)insertions.get(0)).letterDist.length,numFixed)); i++)
		{
			// score codes[] according to the various interactions
			p = 1;
			for(int j = 0; j < interactions.size(); j++)
			{
				p *= ((ClusterInteraction)interactions.get(j)).getSubstProb(codes);
			}
			z += p;
			// push a 1 into codes and carry
			codes[numFixed-1] += 1;
			
			
			int place = numFixed-1;
			while((codes[place] > 3) && (place > 0))// carry
			{
				codes[place] = 0;
				place--;
				codes[place]++;
			}
			
		}
		normalizationZ = z;
		System.out.println("HairpinNode.Normalize normalized a node: " + z);
	}

	public void pretendToNormalize()
	{
		double z = 0;
		double p = 1;
		int[] codes = new int[numFixed];
		
		maxLengths = new int[insertions.size()]; // maximum possible numbers of insertions
		for(int l = 0; l < insertions.size(); l++)
		{
			maxLengths[l] = ((InsertionDistribution)insertions.get(l)).lengthDist.length-1;
		}

		
		for (int i = 0; i < numFixed; i++)
			codes[i] = 0;
				
		normalizationZ = 1;
//		System.out.println("Normalization: "+z);
	}

	/**
	 * This method generates interacting bases for the HairpinNode
	 * @return
	 */
	public String generateInteractingBases()
	{
		double z = 0;
		double p = 1;

		double limit = Math.random() * normalizationZ;
		
		int[] codes = new int[numFixed];
		
		for (int i = 0; i < numFixed; i++)
			codes[i] = 0;

		codes[numFixed - 1] -= 1;
		
		int i = 0;
		
		while((z <= limit) && (i < Math.pow(((InsertionDistribution)insertions.get(1)).letterDist.length,numFixed)))
		{
			// push a 1 into codes and carry
			codes[numFixed-1] += 1;

			int place = numFixed-1;
			while((codes[place] > ((InsertionDistribution)insertions.get(1)).letterDist.length) && (place > 0))// carry
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
//				System.out.println(((ClusterInteraction)interactions.get(j)).getSubstProb(codes));
			}
			z += p;
			
			i++;
			
//			System.out.println(i+"  "+codes[0]+" "+codes[1]+" "+codes[2]+" "+codes[3]+" "+z+" "+limit);
		}
		
//		char letters[] = new char[]{'a','c','g','u'};
		char letters[] = new char[]{'A','C','G','U'};
		
		
// 		take first numLeftInt and return separately from numRightInt
		String leftStr = "";

//		System.out.println(codes[0]+" "+codes[1]+" "+codes[2]+" "+codes[3]+" ");
		
		for(int f = 0; f < numFixed; f++)
			leftStr += letters[codes[f]];
		
		return leftStr;
	}
	
	/**
	 * This method generates insertions for the HairpinNode
	 * @return
	 */
	public Vector generateInsertions()
	{
		Vector ins = new Vector(numFixed);
		String temp = "";
		for(int i = 0; i < ins.capacity(); i++)
		{
			 temp = ((InsertionDistribution)insertions.get(i)).generate();
			 ins.add(temp);
		}
		return ins;
	}
	
	public String generate()
	{
		String bases = generateInteractingBases();
		
		Vector subseq = generateInsertions();
		
		String left = "";
		for(int i = 0; i < numFixed; i++)
			left += bases.charAt(i)+(String)subseq.get(i);

		System.out.println("Hairpin "+left);
		
		return "<" + left + ">";
	}	
	
	
	public void computeMaxLogProb(Sequence seq, int i, int j)
	{
		/**
		 * If the subsequence from i to j is in the range to be parsed by this node, 
		 */
		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax) && (i <= j))
		{
			double p;	    // maximum log probability found so far
			double pnew;    // probability of current insertion possibility
			int a = 0; 		// number of insertions on the left
			int[] intCodes = new int[numFixed];            // codes of interacting bases
			int[] insLengths = new int[insertions.size()]; // one combination of insertion numbers 
			int[] optInsLengths = new int[insertions.size()]; // optimal combination of insertions
			
			for(int l = 0; l < insertions.size(); l++)
			{
                insLengths[l] = 0;
			}
			/*
			 * loop through insertion combinations.
			 * */
			p = -1d/0d;               // start with negative infinity
			p = -1000;                // start with a small probility
			
			while(insLengths[0] <= maxLengths[0])
			{
				// if the number of insertions is correct to fill in the space ...
				if (i+numFixed+a-1 == j)
				{
					// pull out the numFixed interacting base codes
	
					pnew = 0;
	
					int k = i;
					for(int m = 0; m < numFixed; m++)
					{
						intCodes[m] = seq.code[k];
						int[] insCodes = new int[insLengths[m]];
						for (int q = 0; q < insLengths[m]; q++)
							insCodes[q] = seq.code[k+1+q];
						pnew += ((InsertionDistribution)insertions.get(m)).computeLogProb(insCodes);
						k += insLengths[m] + 1;
					}
						
					// score codes[] according to the various interactions
					
					for(int m = 0; m < interactions.size(); m++)
					{
						pnew += ((ClusterInteraction)interactions.get(m)).getLogSubstProb(intCodes);
					}
					
					if(pnew > p)
					{
						//should we keep the codes array for insertions?
						p = pnew;
						for (int z = 0; z < insLengths.length; z++)
						{
							optInsLengths[z] = insLengths[z];
						}
					}
				}
				
				// push a 1 into insLengths and carry
				insLengths[numFixed-1]++;

				int place = numFixed-1;
				while((insLengths[place] > maxLengths[place]) && (place > 0))// carry
				{
					insLengths[place] = 0;
					place--;
					insLengths[place]++;
				}

				// update the total number of insertions
				
				a = 0;
				for(int m = 0; m < numFixed; m++)
				{
					a += insLengths[m];
				}
			
			}	// for loop that sets maxLogProb[i-super.iMin][j-super.jMin]
			
  			maxLogProb[i-iMin][j-jMin] = p;
  			myGen[i-iMin][j-jMin] = new genData(optInsLengths);
		} // end if loop
	}// end method computMaxLogProb

	
	
	public void traceback(int i, int j)
  	{
  		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j))
  		{
  			setOptimalAndRelease(i,j);
  		}
  		else
  		{
  			System.out.println("Hairpin node out of range");
  			System.out.println("Indices "+rightIndex+" "+leftIndex);
  		}
  	}

	
	public String showParse(String n)
  	{
		int[] insert = optimalGen1.splitpoints;
		int i = optimalGen1.i;
		
		String left = n.substring(i,i+1);

		for(int l = 0; l < numFixed-1; l++)
			{
				left += n.substring(i+1,i+1+insert[l]);
				for(int f = 0; f < maxLengths[l] - insert[l] ; f++)
				{
					left += "-";
				}
				i += insert[l]+1;
System.out.println("HairpinNode.showParse "+i+" "+n.length()+" "+numFixed);
                left += n.substring(i,i+1);
			}
		
		return "<" + left + ">";
  	}
	
	
	public String header()
	{
		String left = "";

		for(int l = 0; l < numFixed; l++)
			{
				left = left + "F";
				for(int f = 0; f < maxLengths[l]; f++)
				{
					left += "-";
				}
			}
		
		return "<" + left + ">";
	}
}
