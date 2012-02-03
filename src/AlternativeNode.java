import java.util.LinkedList;

/**
 * This class is used to simulate alternative nodes
 * It works really well.
 * @author meg pirrung
 *
 */
public class AlternativeNode extends BranchingNode {

    int branches;
    double[] probDist;
    double[] logProbDist;
    LinkedList eNodes;
    
	/**
	 * This is the AlternativeNode construtor
	 * @param prev is a pointer to the previous node
	 * @param b is the number of branches for this alternative node
	 * @param p is the array of probabilities to determine which branch should be chosen
	 */
	AlternativeNode(Node prev, int b, double[] p, int lI, int rI)
	{
		super(prev, "AlternativeNode", b, new LinkedList(), lI, rI);
		branches = b;
		probDist = rnaMath.normalize(p);
		logProbDist = rnaMath.log(p);
		eNodes = new LinkedList();
	}
	
	
	public String getParams()
	{
		return super.children.size() + " branches";
	}
	
	
	public String generate()
  	{        
		return ((Node)children.get(rnaMath.rando(probDist))).generate(true);
  	}

	
	void computeMaxLogProb(Sequence seq, int i, int j)
  	{
		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j)) {
			double p = -1d/0d;
			double pnew;
			int branch = 0;
			for(int l = 0; l < probDist.length; l++)
			{
			//  System.out.println("Alternative node: i="+i+" j="+j+" l="+l);
				pnew = logProbDist[l] + ((Node)super.children.get(l)).getMaxLogProb(i,j); 
				if(pnew > p)
				{
					p = pnew;
					branch = l;
				}
			}
			maxLogProb[i-iMin][j-jMin] = p;
			myGen[i-iMin][j-jMin] = new genData(branch);			
		}
	}
	
	
  	
	public void traceback(int i, int j)
  	{
  		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j))
  		{
  			setOptimalAndRelease(i,j);

  			for (int k = 0; k < branches; k++)
  	   			((Node)super.children.get(k)).traceback(i,j);
  				
  			
//   			int branch = optimalGen.branch;
//   			((Node)super.children.get(branch)).traceback(i,j);
  		}
  		else
  			System.out.println("Alternative node out of range");
  	}
   	
  	
	public String showParse(String n)
  	{
  			String parse = "";
  			for (int i = 0; i < branches; i++)
  			{
  				String newChild = ((Node)children.get(i)).showParse(n);
  				if (i != optimalGen1.branch)
  				{
//  				newChild = newChild.replace("-","+");
  					newChild = newChild.replace("A","+");
  					newChild = newChild.replace("C","+");
  					newChild = newChild.replace("G","+");
  					newChild = newChild.replace("U","+");
  					newChild = newChild.replace("a","+");
  					newChild = newChild.replace("c","+");
  					newChild = newChild.replace("g","+");
  					newChild = newChild.replace("u","+");
  				}
  				parse = parse + i + newChild + i;
  			}
  			
			return parse;
  	}
  	
  	
	public String header()
  	{		
			String parse = "";
  			for (int i = 0; i < branches; i++)
  			{
  				String newChild = ((Node)children.get(i)).header();
  				parse = parse + i + newChild + i;
  			}
  			
			return parse;
  	}

}
