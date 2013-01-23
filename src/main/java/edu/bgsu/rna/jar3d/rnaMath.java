package edu.bgsu.rna.jar3d;
/**
 * @author meg pirrung
 *
 */
public class rnaMath {

	/**
	 * This class contains math functions for use specifically in the Parser
	 * project
	 */
	public rnaMath() {
		super();
	}
	
	/**
	 * This method selects an array index based on probabilities included in the arrray
	 * @param dist is a distribution matrix
	 * @return i
	 */
	public static int rando(double[] dist)
	{
		int n = dist.length;
		double u = Math.random();
		int i = 0;
		double s = dist[0];
		
		while((u>s)&&(i<n-1))
		{
			i++;
			s = s+dist[i];
		}
		return i;
	}
	
	/**
	 * This method selects an array index based on probabilities included in the arrray
	 * @param dist is a distribution matrix
	 * @return i
	 */
	public static int[] randod(double[][] dist)
	{
		int n = dist[0].length;
		int m = dist.length;
		double u = Math.random();
		int i = 0;
		int j = 0;
		double s = dist[0][0];
		
		while((u>s)&&(i<n-1)&&(j<m-1))
		{
			i++;
			j++;
			s = s+dist[j][i];
		}
		return new int[]{j,i};
	}
	
	/**
	 * This method checks to see if all array entries add up to 1, if not it fixes that
	 * @param dist is a distribution matrix
	 * @return
	 * NOTE : Normalization has been disabled!
	 */
	public static double[] normalize(double[] dist)
	{
		double newdist[] = new double[dist.length];
		double sum = 0;
		for(int i = 0; i < dist.length; i++)
		{
			sum += dist[i];
		}
		for(int i = 0; i < dist.length; i++)
		{
			//Does NOT normalize!
			newdist[i] = dist[i];
		}
		return newdist;
	}
	
	/**
	 * This method checks to see if all array entries add up to 1, if not it fixes that
	 * @param dist is  a doubly indexed distribution matrix
	 * @return
	 * NOTE : Normalization has been disabled!
	 */
	public static double[][] normalize(double[][] dist)
	{
		double newdist[][] = new double[dist.length][dist.length];
		double sum = 0;
		for(int i = 0; i < dist.length; i++)
		{
			for(int j = 0; j < dist.length; j++)
			{
				sum += dist[i][j];
			}
		}
		for(int i = 0; i < dist.length; i++)
		{
			for(int j = 0; j < dist.length; j++)
			{
				//Does NOT normalize!
				newdist[i][j] = dist[i][j];
			}
		}
		return newdist;
	}

	/**
	 * This method takes the logarithm of an array of doubles
	 * @param ary
	 * @return
	 */
	public static double[] log(double[] ary)
	{
		double newary[] = new double[ary.length];
		
		for(int i = 0; i < ary.length; i++)
		{
			newary[i] = Math.log(ary[i]);
		}
		return newary;
	}
	
	/**
	 * This method takes the logarithm of a doubly indexed array
	 * of doubles
	 * @param ary
	 * @return
	 */
	public static double[][] logd(double[][] ary)
	{
		double newary[][] = new double[ary.length][ary.length];
		
		for(int i = 0; i < ary.length; i++)
		{
			for(int j = 0; j < ary.length; j++)
			{
				newary[i][j] = Math.log(ary[i][j]);
			}
		}
		return newary;
	}
}
