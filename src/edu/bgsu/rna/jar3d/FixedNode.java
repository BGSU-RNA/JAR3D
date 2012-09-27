package edu.bgsu.rna.jar3d;
/**
 * This class is used to code one single fixed nucleotide
 * @author Craig Zirbel
 *
 */
public class FixedNode extends BasicNode {
	InsertionDistribution lInsDist, rInsDist;
	
	// arrays of logs of all probabilities

    int alt;
    int interacttype;
    int Position;                          // 1 is left, 2 is right
   
	/**
	 * This is the InitialNode constructor
	 * @param prev contains a pointer to the previous node
	 * @param lLenDist is the left length distribution matrix
	 * @param lLetDist is the left letter distribution matrix
	 * @param rLenDist is the right length distribution matrix
	 * @param rLetDist is the right letter distribution matrix
	 */
	FixedNode(Node prev, double dProb, double[] letDist, int posi, int lI, int rI)
	{
		super(prev, "FixedNode", lI, rI);

		Position = posi;
		
		double probs[] = new double[2];
		probs[0] = dProb;
		probs[1] = 1-dProb;
		
		double zeroLength[] = new double[1];
		zeroLength[0] = 1;
		
		if (Position == 1) {
			lInsDist = new InsertionDistribution(probs, letDist);
			rInsDist = new InsertionDistribution(zeroLength, letDist);
		}
		else {
			lInsDist = new InsertionDistribution(zeroLength, letDist);
			rInsDist = new InsertionDistribution(probs, letDist);
		}
	}

  	public String generate()
  	{       
  		return lInsDist.generate() + super.child.generate(true) + rInsDist.generate();
  	}

	
	public void computeMaxLogProb(Sequence seq, int i, int j)
	{
		/**
		 * If the subsequence from i to j is in the range to be parsed by this node, 
		 */
		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax) && (i <= j))
		{
			double p;	// maximum log probability found so far
			double pll;	// contribution to total log prob from left length
			double pli;	// contribution to total log prob from left inserted letters
			double prl; // contribution to total log prob from right length
			double pri; // contribution to total log prob from right inserted letters
			double priarray[] = new double[rInsDist.logLengthDist.length];
			double pnew; // probability of current insertion possibility
			int a; 		// number of insertions on the left
			int b; 		// number of insertions on the right
			int aa = 0, bb = 0;

			
			p = -1d/0d;               // start with negative infinity
			// for loop that sets maxLogProb[i-super.iMin][j-super.jMin]
			pli = 0;					// 0 left insertions so far

			priarray[0] = 0;
			for(b = 1; b < Math.min(j-i+1,rInsDist.lengthDist.length); b++)
  			{
				priarray[b] = priarray[b-1];
				int c = seq.code[j-b+1];
				priarray[b] += rInsDist.logLetterDist[c];
  			} // end for loop

			for(a = 0; a < Math.min(j-i+1,lInsDist.logLengthDist.length); a++)
  	  		{
  				pll = lInsDist.logLengthDist[a];
  				if (a > 0) {
  					pli = pli + lInsDist.logLetterDist[seq.code[i+a-1]];
  				}
  				pri = 0;				// 0 right insertions so far
  				for(b = 0; b < Math.min(j-i+1-a,rInsDist.lengthDist.length); b++)
  	  			{
  					prl = rInsDist.logLengthDist[b];
  					pri = priarray[b];
					pnew = super.child.getMaxLogProb(i+a,j-b);
  	  				pnew = pnew + pll + pli + prl + pri;
  	  				if (pnew > p) {
  	  					p = pnew;
  	  					aa = a;
  	  					bb = b;
  	  				} // end if
  	  			} // end for loop
  	  		} // end for loop
  			maxLogProb[i-iMin][j-jMin] = p;
  			myGen[i-iMin][j-jMin] = new genData(aa,bb);
		} // end if loop
	}// end method computMaxLogProb

	
  	
	public void traceback(int i, int j)
  	{
  		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j))
  		{
   			setOptimalAndRelease(i,j);
  			int a = optimalGen1.numLeftIns;
  			int b = optimalGen1.numRightIns;
  			super.child.traceback(i+a,j-b);
  		}
  		else
  		{
  			System.out.println("FixedNode.traceback: Fixed node out of range");
  			System.out.println("  i="+i+" j="+j+" iMin="+super.iMin+" iMax="+super.iMax+" jMin="+super.jMin+" jMax="+super.jMax);
  			System.out.println("  Indices "+leftIndex+" "+rightIndex);
  			System.out.println("  Previous was " + super.previous.mytype + " iMax=" + super.previous.iMax + " jMin=" + super.previous.jMin);
  		}
  	}
  	
	public String showParse(String n)
  	{
  			int a = optimalGen1.numLeftIns;
  			int b = optimalGen1.numRightIns;
  			int i = optimalGen1.i;
  			int j = optimalGen1.j;
  			
  			String left = n.substring(i,i+a);
  			int lSize = left.length();
  			for(int f = 0; f < lInsDist.lengthDist.length - lSize; f++)
  				left += "-";
  			
  			String right = n.substring(j-b+1,j+1);
  			int rSize = right.length();
  			for(int g = 0; g < rInsDist.lengthDist.length - rSize; g++)
  				right = "-"+right;
  			
  		    return left + super.child.showParse(n) + right;
  	}
  	
  	
	public String header()
  	{		
  			String left = "F";
  			for(int f = 0; f < lInsDist.lengthDist.length-1; f++)
  				left += "-";
  			
  			String right = "F";
  			for(int g = 0; g < rInsDist.lengthDist.length-1; g++)
  				right = "-"+right;
  			
  		    return left + super.child.header() + right;
  	}

	public String showCorrespondences(String letters)
  	{
  			int a = optimalGen1.numLeftIns;
  			int b = optimalGen1.numRightIns;
  			int i = optimalGen1.i;
  			int j = optimalGen1.j;

  			String left = "";
  			for(int k = i; k < i + a; k++)
  				left += "SSS_Position_" + (k+1) + "_" + letters.charAt(k) + " JAR3D_aligns_to " + "MMM_Node_" + number + "_Position_1" + "\n";
  			
  			String right = "";
  			for(int k = j-b+1; k <= j; k++)
  				right = "SSS_Position_" + (k+1) + "_" + letters.charAt(k) + " JAR3D_aligns_to " + "MMM_Node_" + number + "_Position_2" + "\n" + right;
 			
  			return left + super.child.showCorrespondences(letters) + right;
  	}
} // end class FixedNode