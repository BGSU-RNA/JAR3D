package edu.bgsu.rna.jar3d;
import java.util.*;

/**
 * This class is used to simulate Junctions
 * @author meg pirrung
 *
 */
public class JunctionNodeMultipleJunctions extends BranchingNode {

    int branches;
//    long start, stop, elapsed;
    
	/**
	 * This is the JunctionNode constructor
	 * @param prev contains a pointer to the previous node
	 * @param b is the number of branches in the junction
	 */
	JunctionNodeMultipleJunctions(Node prev, int b, int lI, int rI)
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

	
	void computeMaxLogProbMultiple(Sequence seq, int i, int j)
	{
//		start = System.currentTimeMillis();
		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j)) 
		{		
			int[] splitPoints = new int[branches+1];
			int[] minSplit = new int[branches+1];
			int[] maxSplit = new int[branches+1];
			int[] optSplitPoints = new int[branches+1];
			double p;	// maximum log probability found so far
			double pnew; // probability of current insertion possibility
						
			minSplit[0] = i-1;
			maxSplit[0] = i-1;
			
			minSplit[branches] = j;
			maxSplit[branches] = j;
			
			for(int f = 1; f < branches; f++)
			{
				minSplit[f] = ((Node)super.children.get(f)).iMin-1;
				maxSplit[f] = ((Node)super.children.get(f)).iMax-1;
			}
			
			for(int f = 1; f < branches; f++)
			{
				minSplit[f] = Math.max(minSplit[f], ((Node)super.children.get(f-1)).jMin + 1);
				maxSplit[f] = Math.min(maxSplit[f], ((Node)super.children.get(f-1)).jMax + 1);
			}
			
			for(int f = 0; f < branches+1; f++)
			{
				splitPoints[f] = minSplit[f];
				optSplitPoints[f] = minSplit[f];
			}
			
			
			p = -1d/0d;               // start with negative infinity

			while(splitPoints[0] <= maxSplit[0])
			{
				pnew = 0;
				for(int n = 0; n < branches; n++)
					pnew += ((Node)super.children.get(n)).getMaxLogProb(splitPoints[n]+1, splitPoints[n+1]);
			
				if(pnew > p)
				{
					p = pnew;
					for (int z = 0; z < branches+1; z++)
					{
						optSplitPoints[z] = splitPoints[z];
					}
				}
					
				// push a 1 into insLengths and carry
				splitPoints[branches-1]++;

				int place = branches-1;
				while((splitPoints[place] > maxSplit[place]) && (place > 0))// carry
				{
					splitPoints[place] = Math.max(minSplit[place], splitPoints[place-1]+1);
					place--;
					splitPoints[place]++;
				}

			}
			// for loop that sets maxLogProb[i-super.iMin][j-super.jMin]
			
  			maxLogProb[i-iMin][j-jMin] = p;
  			myGen[i-iMin][j-jMin] = new genData(optSplitPoints);
		} // end if loop
//		stop = System.currentTimeMillis();
//		elapsed = stop - start;
		//System.out.println("Junction node - Compute max log prob time: " + elapsed+ " milliseconds");
	}
	
	
  	
	public void traceback(int i, int j)
  	{
  		if ((i >= super.iMin) && (i <= super.iMax) && (j >= super.jMin) && (j <= super.jMax)  && (i <= j))
  		{
  			setOptimalAndRelease(i,j);

			((Node)super.children.get(0)).traceback(optimalGen1.splitpoints[0], optimalGen1.splitpoints[1]);
			for(int n = 1; n < branches; n++)
				((Node)super.children.get(n)).traceback(optimalGen1.splitpoints[n]+1, optimalGen1.splitpoints[n+1]);
	}
  		else
  		{
  			System.out.println("Junction node out of range");
  			System.out.println("Indices "+leftIndex+" "+rightIndex);
  		}
  	}
   	
  	
	public String showParse(String n)
  	{
		String parse = "";
		for(int f = 0; f < branches; f++)
			parse += ((Node)super.children.get(f)).showParse(n);
		return parse;
  	}
  	
  	
	public String header()
  	{
  			String header = "";
  			for(int f = 0; f < branches; f++)
  				header += ((Node)super.children.get(f)).header();
			return header;
  	}
}// end class JunctionNode