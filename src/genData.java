
/**
 * This method contains the maximum proability way that each node generated a certain
 * sequence of letters
 * @author meg pirrung
 */
public class genData {
	int numLeftIns, numRightIns;
	char lPair, rPair;
	boolean deleted;
	String genType;
    int splitpoints[];
    int insLengths[];
    public int i;            // maybe this should be a Node variable?
    public int j;
    public int branch;
    
    /**
     * This is the genData for InitialNodes
     * @param a is the number of insertions on the left
     * @param b is the number of insertions on the right
     */
    genData(int a, int b)
	{
		numLeftIns = a;
		numRightIns = b;
	}
    
    /**
     * This is the genData for BasepairNodes
     * @param m is the left base
     * @param n is the right base
     * @param a is the number of insertions on the left
     * @param b is the number of insertions on the right
     */
    genData(char m, char n, int a, int b)
    {
    	lPair = m;
    	rPair = n;
    	numLeftIns = a;
    	numRightIns = b;
    	deleted = false;
    }
	
    /**
     * This is the genData for deleted BasepairNodes and ClusterNodes
     * @param del
     */
    genData(boolean del)
    {
    	deleted = true;
    }
    
    /**
     * This is the genData for HairpinNodes
     * @param type is the type of genData
     */
    genData(String type)
	{
		genType = type;
	}

    /**
     * This is the genData for JunctionNodes with more than one branch
     * @param a
     */
    genData(int[] a)
	{
    	splitpoints = a;
	}
    
    /**
     * This is the genData for ClusterNodes
     * @param L is a matrix of the length of insertions for each position
     * @param aa is the number of left insertions
     * @param bb is the number of right insertions
     */
    genData(int[] L,int aa, int bb)
	{
   		insLengths = L;
   		numLeftIns = aa;
   		numRightIns = bb;
	}

    /**
     * This is the genData for AlternativeNodes
     * @param b is the branch to be used
     */
    genData(int b)
    {
    	branch = b;
    }
    
}
