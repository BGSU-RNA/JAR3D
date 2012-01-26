/**
 * Node is the wrapper class for all other node classes
 * @author meg pirrung
 *
 */
class Node {
	public String mytype;
	public Node next, previous, child;
	public int iMin, iMax, jMin, jMax;
	public int maxlength;
	public double maxLogProb[][];
	public genData myGen[][];
	public double optimalMaxLogProb;
	public genData optimalGen1, optimalGen2;
	public boolean visible, active;
	public int rightIndex, leftIndex;
	public double currentMaxLogProb;
	
	/**
	 * This is the Node constructor
	 * @param prev contains a pointer to the previous node
	 * @param type contains the type of node
	 */
	Node(Node prev, String type, int lI, int rI)
	{
		mytype = type;
		next = null;
		child = null;
		previous = prev;
		visible = true;
		active = false;
		rightIndex = rI-1;            // adjust between Matlab convention and Java convention
		leftIndex = lI-1;
	}
	
	Node()
	{
		mytype = null;
		next = null;
		child = null;
		previous = null;
		visible = true;
		active = false;
		rightIndex = 0;
		leftIndex = 0;
	}
	
	/**
	 * This method returns the type of node
	 * @return
	 */
	public String getType()
	{
		return mytype;
	}

	/**
	 * This is the overloaded generate method that is rewritten differently
	 * in each node
	 * @param del tells whether this node is deleted or not
	 * @return A string of what this node generated
	 */
	public String generate(boolean del)
	{
		return "AG";
	}
	
	/**
	 * This method sets the minimum and maximum number of probabilities stored
	 * @param a 
	 * @param b
	 * @param c
	 * @param d
	 */
	public void setMinMax(int range, Sequence seq)
	{		
		/*
		System.out.println("Node.setMinMax");
		System.out.println(this.mytype);
		System.out.println(seq.cti.length);
		System.out.println(seq.itc.length);
		System.out.println(seq.ctiFirst.length);
		System.out.println(seq.itcFirst.length);
		
		
		for(int c = 0; c < seq.letters.length(); c++)
		{
			System.out.print(seq.letters.charAt(c));
		}

		System.out.println();

		System.out.print("cti: ");
		for(int c = 0; c < Math.min(20,seq.cti.length); c++)
		{
			System.out.print(seq.cti[c]+" ");
		}

		System.out.println();
		
		System.out.print("itc: ");
		for(int c = 0; c < Math.min(20,seq.itc.length); c++)
		{
			System.out.print(seq.itc[c]+" ");
		}

		System.out.println();
		*/
		
		int leftCol = seq.itcFirst[Math.min(leftIndex,seq.itcFirst.length-1)];    // 
		leftCol = Math.min(leftCol, seq.cti.length-1);
		//int leftCol = seq.itcFirst[leftIndex];    // 
		iMin = seq.cti[leftCol] - range;
		iMax = seq.cti[leftCol] + range;
		int rightCol = seq.itcFirst[Math.min(rightIndex,seq.itcFirst.length-1)];
		rightCol = Math.min(rightCol, seq.cti.length-1);
		jMin = seq.cti[rightCol] - range;
		jMax = seq.cti[rightCol] + range;
		
		iMin = Math.max(iMin,0);
		iMax = Math.min(iMax,seq.nucleotides.length());
		jMin = Math.max(jMin,0);
		jMax = Math.min(jMax,seq.nucleotides.length());

		if((iMax > iMin) && (jMax > jMin))
		{
			maxLogProb = new double[iMax-iMin+1][jMax-jMin+1];
			myGen = new genData[iMax-iMin+1][jMax-jMin+1];
		}
		else
		{//look here
			System.out.println("This sequence is too short for the model and the limits on iMin and iMax");
			System.out.println(leftIndex+" "+rightIndex);
			System.out.println(iMin+" "+iMax+" "+jMin+" "+jMax+" ");
			System.out.println(mytype);
			maxLogProb = new double[1][1];
			myGen = new genData[1][1];
		}
	}
	
	/**
	 * This is the overloaded method that computes the logarithm of the
	 * maxium probability way that a node generated a certain sequence
	 * @param seq
	 * @param i
	 * @param j
	 */
	void computeMaxLogProb(Sequence seq, int i, int j)
	{}
	
	/**
	 * This method returns the logarithm of the maximum probability
	 * @param i
	 * @param j
	 * @return
	 */
	public double getMaxLogProb(int i, int j)
	{
		if ((i >= iMin) && (i <= iMax) && (j >= jMin) && (j <= jMax))
		{
			if (maxLogProb[i-iMin][j-jMin] > currentMaxLogProb)
				currentMaxLogProb = maxLogProb[i-iMin][j-jMin];
			return maxLogProb[i-iMin][j-jMin];
		}
		else
		{
			return -1d/0d; // return negative infinity
		}
	}

	String getParams()
	{return"";}

	public void InternalCalculations()
	{}
	
	/**
	 * This method is overloaded for each different node class
	 * @param i
	 * @param j
	 */
	public void traceback(int i, int j)
	{}

	/**
	 * This method is overloaded for each different node class
	 * It scores a string of characters n to the model represented by the nodes in
	 * this sequence
	 * @param n
	 * @return
	 */
	public String showParse(String n)
  	{
		return "Blank";
  	}
	
	/**
	 * This method is used for creating an alignment header
	 * @return
	 */
	public String header()
	{
		return " ";
	}

	/**
	 * This method keeps the optimal maxLogProb and releases the
	 * extra array space
	 * @param i
	 * @param j
	 */
	public void setOptimalAndRelease(int i, int j)
	{
		// if this hasn't already been done, then do this:
		optimalMaxLogProb = maxLogProb[i-iMin][j-jMin];
		
		
		if (optimalMaxLogProb < -99999999)
		{
/*			System.out.print("Node.setOptimalAndRelease "+mytype+" "+iMin+"<="+i+"<="+iMax+" "+jMin+"<="+j+"<="+jMax);
			System.out.print(" length "+myGen.length);
			System.out.print(" leftIndex "+leftIndex);
			System.out.print(" rightIndex "+rightIndex);
			System.out.println(" maxlogprob " + optimalMaxLogProb);
*/		}
		
		optimalGen1 = myGen[i-iMin][j-jMin];
		optimalGen1.i = i;
        optimalGen1.j = j;
		myGen = new genData[0][0];
		maxLogProb = new double[0][0];
	}
}
