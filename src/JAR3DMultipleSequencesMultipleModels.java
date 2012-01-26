
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;

/**
 * This is a driver program
 * @author meg pirrung
 *
 */
public class JAR3DMultipleSequencesMultipleModels {
	public static void main(String[] args) {

		int numSequences = 1000;                            // make sure this is larger than needed

		String loopType = "HL";
		Vector modnames = Sequence.getModelNames(loopType);
        String FASTAName = "";
		Vector sData;
		Vector scores = new Vector();
		
        for (int m=0; m<modnames.size();m++)                // loop through sets of sequences
        {
            FASTAName = ((String)modnames.get(m)).replace(".txt",".fasta");
            System.out.println("Aligning sequences from "+FASTAName);
//            System.out.println("JAR3DMultipleSequencesMultipleModels: sequences\\"+FASTAName);
            sData = Alignment.loadFasta(FASTAName); 
        	String newscores = Alignment.getSortedHLAlignment(sData,modnames,numSequences,100);
            scores.add(newscores);
        }
        
        System.out.println("");
        System.out.print("Names = {");
        for (int m=0; m < modnames.size();m++)
        {
        	System.out.print("'"+((String)modnames.get(m)).replace(".txt","")+"'");
        	if (loopType.equals("IL")) 
        		System.out.print(",'R "+((String)modnames.get(m)).replace(".txt","")+"'");
            if (m < (modnames.size()-1))
            	System.out.print(",");
        }
        System.out.println("};");
        
        for (int m=0; m < modnames.size();m++)
        {
        	System.out.println("S("+(m+1)+",:) = ["+(String)scores.get(m)+"];");
        }
}
	
}
