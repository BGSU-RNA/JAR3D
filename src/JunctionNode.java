import java.util.*;

/**
 * This class is used to simulate Junctions
 * @author meg pirrung
 *
 */
public class JunctionNode extends BranchingNode {

    int branches;
    long start, stop, elapsed;
    
	/**
	 * This is the JunctionNode constructor
	 * @param prev contains a pointer to the previous node
	 * @param b is the number of branches in the junction
	 */
	JunctionNode(Node prev, int b, int lI, int rI)
	{
		super(prev, "JunctionNode", b, new LinkedList(), lI, rI);
		branches = b;
	}
	
	
	public String getParams()
	{
		return super.children.size() + " branches";
	}
	
	
	public String generate()
  	{        
		String bGen = "";
		for(int i = 0; i < children.size(); i++)
			bGen += ((Node)children.get(i)).generate(true); // true is a placeholder
		return bGen;
  	}

	
	void computeMaxLogProb(Sequence seq, int i, int j)
  	{
		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j)) {
			double p = -1d/0d;
			double pnew;
			int kk = 0, mm = 0;
			switch(branches)
			{
				case 2:
					for(int k = i; k < j; k++)
					{		   // branch one, from i to k        branch two, from k+1 to j
							   // i-------------k--------------j
						pnew = ((Node)super.children.get(0)).getMaxLogProb(i,k) + ((Node)super.children.get(1)).getMaxLogProb(k+1,j);
						if (pnew > p)
						{
							p = pnew;
							kk = k;
						}
					}				
					maxLogProb[i-iMin][j-jMin] = p;
					myGen[i-iMin][j-jMin] = new genData(new int[]{kk});
					/*
					if ((p < -99999999))
					{
						System.out.println("Junction node with -Inf prob.  "+i+" "+j+" "+leftIndex+" "+rightIndex);
					}
					*/
				break;
			case 3:
				for(int k = i; k < j; k++)
				{
					for(int m = k+1; m < j; m++)
					{
						       // branch one from i to k         branch two from k+1 to m          branch three from m+1 to j
							   // i-------------k----m---------j
						pnew = ((Node)super.children.get(0)).getMaxLogProb(i,k) + ((Node)super.children.get(1)).getMaxLogProb(k+1,m) + ((Node)super.children.get(1)).getMaxLogProb(m+1,j);
						if (pnew > p)
						{
							p = pnew;
							kk = k;
							mm = m;
						}
					}
				}		
				maxLogProb[i-iMin][j-jMin] = p;
				myGen[i-iMin][j-jMin] = new genData(new int[]{kk,mm});
				break;
			}
		}
	}
	
	
  	
	public void traceback(int i, int j)
  	{
  		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j))
  		{
  			setOptimalAndRelease(i,j);

   			int kk = optimalGen1.splitpoints[0];
			switch(branches)
			{
				case 2:
					((Node)children.get(0)).traceback(i,kk);
					((Node)children.get(1)).traceback(kk+1,j);
					break;
				case 3:
					int mm = optimalGen1.splitpoints[1];
					((Node)children.get(0)).traceback(i,kk);
					((Node)children.get(1)).traceback(kk+1,mm);
					((Node)children.get(2)).traceback(mm+1,j);
					break;
			}
  		}
  		else
  			System.out.println("Out of range");
  	}

/*	
	public String showParse(String n, int i, int j)
  	{
   		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j))
  		{
   			String parse = "";
   			setOptimalAndRelease(i,j);
   			int kk = optimalGen.splitpoints[0];
			switch(branches)
			{
				case 2:
					parse = ((Node)children.get(0)).showParse(n,i,kk);
					parse += ((Node)children.get(1)).showParse(n,kk+1,j);
					break;
				case 3:
					int mm = optimalGen.splitpoints[1];
					parse = ((Node)children.get(0)).showParse(n,i,kk);
					parse += ((Node)children.get(1)).showParse(n,kk+1,mm);
					parse += ((Node)children.get(2)).showParse(n,mm+1,j);
					break;
			}
			return parse;
  		}
   		else
   			return "Out of range";
   		
  	}
*/	
   	
  	
	public String showParse(String n)
  	{
  			String parse = "";
			switch(branches)
			{
				case 2:
					parse = ((Node)children.get(0)).showParse(n);
					parse += ((Node)children.get(1)).showParse(n);
					break;
				case 3:
					parse = ((Node)children.get(0)).showParse(n);
					parse += ((Node)children.get(1)).showParse(n);
					parse += ((Node)children.get(2)).showParse(n);
					break;
			}
			return parse;
  	}
  	
  	
	public String header()
  	{
  			String parse = "";
			switch(branches)
			{
				case 2:
					parse = ((Node)children.get(0)).header();
					parse += ((Node)children.get(1)).header();
					break;
				case 3:
					parse = ((Node)children.get(0)).header();
					parse += ((Node)children.get(1)).header();
					parse += ((Node)children.get(2)).header();
					break;
			}
			return parse;
  	}

}// end class JunctionNode