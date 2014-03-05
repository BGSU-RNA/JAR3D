package edu.bgsu.rna.jar3d;

/**
 * This class holds common parameters and methods used in ClusterNodes
 * @author meg pirrung
 *
 */
public class ClusterInteraction {
	int base1, base2;
	double[][] subsProb;
	double[][] logSubsProb;
	
	/**
	 * This is the ClusterInteraction constructor
	 * @param a is the position of the left side base within the cluster
	 * @param b is the position of the right side base within the cluster
	 * @param subs is a matrix containing substitution probabilities
	 */
	public ClusterInteraction(int a, int b, double[][]subs)
	{
		if (a >= 0 & b >= 0)
		{
			base1 = a;
			base2 = b;
			subsProb = rnaMath.normalize(subs);
			logSubsProb = rnaMath.logd(subsProb);
		}
		else
			System.out.println("ClusterInteraction: Base position is negative" + a + " " + b);
	}
	
	/**
	 * This method returns the substitution probability of the given bases
	 * @param baseString contains the two bases
	 * @return
	 */
	public double getSubstProb(int[] baseString)
	{
        int code1 = baseString[base1];
        int code2 = baseString[base2];
		return subsProb[code1][code2];
	}
	
	/**
	 * This method returns the logarithm of the substitution probability of the given bases
	 * @param baseString contains the two bases
	 * @return
	 */
	public double getLogSubstProb(int[] baseString)
	{
		int code1 = baseString[base1];
		int code2 = baseString[base2];
		return logSubsProb[code1][code2];
	}

}
