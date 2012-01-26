import java.util.*;

/**
 * This is a driver program
 * @author meg pirrung
 *
 */
public class JAR3DGenerateAndParse {
	public static void main(String[] args) {
		
		int modelNum = 5;
		int numSequences = 3;
		
		Vector sequenceData = Alignment.generateSyntheticData(numSequences,modelNum);
		Vector parseData = Alignment.doParse(sequenceData,modelNum,numSequences, " ");
		Alignment.displayAlignment(parseData,sequenceData);

	}
}
