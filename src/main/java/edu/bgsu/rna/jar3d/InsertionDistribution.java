package edu.bgsu.rna.jar3d;

/**
 * This class is used to manage insertion distributions and their common
 * methods
 * @author meg pirrung
 *
 */
public class InsertionDistribution {
	
	double[] lengthDist;
	double[] letterDist;
	double[] logLengthDist;
	double[] logLetterDist;
	
	/**
	 * This is the InsertionDistribution Constructor
	 * @param lDist is the length distribution matrix
	 * @param letDist is the letter distribution matrix
	 */
	public InsertionDistribution(double[] lDist, double[] letDist)
	{
		lengthDist = rnaMath.normalize(lDist);
		letterDist = rnaMath.normalize(letDist);
		
		logLengthDist = rnaMath.log(lengthDist);
		logLetterDist = rnaMath.log(letterDist);
	}
	
	/**
	 * This method generates an insertion based on the length and letter
	 * distribution matrices
	 * @return
	 */
	public String generate()
	{
		char letters[] = new char[]{'A','C','G','U'};

		int length = rnaMath.rando(lengthDist);
		String insert = "";
		for(int i = 0; i < length; i++)
		{
			insert += letters[rnaMath.rando(letterDist)];
		}
		return insert;
	}	
	
	/**
	 * This method computes the logarithm of the probability of a certain
	 * sequence being generated based on the stored length and letter 
	 * distribution matrices
	 * @param codes
	 * @return
	 */
	double computeLogProb(int[] codes)
	{
		double p;
		p = logLengthDist[codes.length];       // probability of this length of insertion
		for (int i = 0; i < codes.length; i++)
		{
			p += logLetterDist[codes[i]];
		}
			return p;
	}

	/**
	 * This method computes the logarithm of the probability of a certain
	 * sequence being generated based on the stored length and letter 
	 * distribution matrices
	 * @param codes
	 * @return
	 */
	double computeProb(int[] codes)
	{
		double p;
		p = lengthDist[codes.length];       // probability of this length of insertion
		for (int i = 0; i < codes.length; i++)
		{
			p *= letterDist[codes[i]];
		}
			return p;
	}

	/**
	 * This method returns the maximum number of insertions possible
	 * @return
	 */
	int maxNumInsertions()
	{
		return lengthDist.length;
	}
}
