package edu.bgsu.rna.jar3d;

import java.io.*;
import java.util.*; 

public class webJAR3D {
	
	public static String double2String(double number, int decimals){
		String S = Double.toString(number);
		int decloc = S.indexOf(".");
		int currentDecs = S.substring(decloc).length()-1;
		if(currentDecs == decimals) return S;
		else if(currentDecs > decimals) return S.substring(0,decloc+decimals+1);
		else {
			int n = decimals - currentDecs;
			for(int i = 0; i < n; i++) S = S+"0";
			return S;
		}
	}
	
	//getQuantiles for internal use by JAR3DMatlab or WebJAR3D functions
	public static double[] getQuantilesA(double[] Scores, String groupNum, String loopType) {
		int n = Scores.length;
		double[] quantiles = new double[n];
		Vector modDist = new Vector();
		Vector modVals = new Vector();
		
		char fsep = File.separatorChar;
		String curDir = System.getProperty("user.dir");
		String modDistFN = curDir + fsep + "Models" + fsep + "Emp. Distributions" + fsep + loopType + "_" + groupNum + ".txt";
		File distFile = new File(modDistFN);
		
		try {
			FileReader inStream = new FileReader(distFile);
			BufferedReader in = new BufferedReader(inStream);
			String lineS;
			String ValueS;
			String DistS;
			double Value;
			double Dist;
			int BreakP;
			while((lineS = in.readLine()) != null){
				BreakP = lineS.indexOf(" ");
				ValueS = lineS.substring(0, BreakP);
				DistS = lineS.substring(BreakP+1);
				Dist = Double.parseDouble(DistS);
				modDist.add(Dist);
				Value = Double.parseDouble(ValueS);
				modVals.add(Value);
			}
		} catch(Exception e) {
			System.out.println("Error reading file " + modDistFN + "\n"+e);
		}
		
		int DistLength = modVals.size();
		
		for(int i = 0; i < n; i++){
			int found = 0;
			int j = DistLength/2;
			int step = j;
			double current;
			double next;
			while(found==0){
				if(j==DistLength-1) {
					quantiles[i] = 1;
					found = 1;
				}
				else if(j==0){
					current = ((Double)modVals.get(j)).doubleValue();
					if(Scores[i] == current || Scores[i] > current){
						quantiles[i] = ((Double)modDist.get(j)).doubleValue();
						found = 1;
					}else{
						quantiles[i] = 0;
						found = 1;
					}
				}else{
					current = ((Double)modVals.get(j)).doubleValue();
					next = ((Double)modVals.get(j+1)).doubleValue();
					if(Scores[i] == current || (Scores[i] > current && Scores[i] < next)){
						quantiles[i] = ((Double)modDist.get(j)).doubleValue();
						found = 1;
					}else if(Scores[i] <current){
						step = Math.max(step/2,1);
						j = j - step;
					}else {
						step = Math.max(step/2,1);
						j = j + step;
					}
				}
			}
		}
		return quantiles;
	}
	//getQuantiles for external use (to be called from Matlab)
	public static double[] getQuantilesB(double[] Scores, String groupNum, String loopType) {
		int n = Scores.length;
		double[] quantiles = new double[n];
		Vector modDist = new Vector();
		Vector modVals = new Vector();
		
		char fsep = File.separatorChar;
//		String curDir = System.getProperty("user.dir");
		String modDistFN = "Models" + fsep + "Emp. Distributions" + fsep + loopType + "_" + groupNum + ".txt";
		File distFile = new File(modDistFN);
		
		try {
			FileReader inStream = new FileReader(distFile);
			BufferedReader in = new BufferedReader(inStream);
			String lineS;
			String ValueS;
			String DistS;
			double Value;
			double Dist;
			int BreakP;
			while((lineS = in.readLine()) != null){
				BreakP = lineS.indexOf(" ");
				ValueS = lineS.substring(0, BreakP);
				DistS = lineS.substring(BreakP+1);
				Dist = Double.parseDouble(DistS);
				modDist.add(Dist);
				Value = Double.parseDouble(ValueS);
				modVals.add(Value);
			}
		} catch(Exception e) {
			System.out.println("Error reading file " + modDistFN + "\n"+e);
		}
		
		int DistLength = modVals.size();
		
		for(int i = 0; i < n; i++){
			int found = 0;
			int j = DistLength/2;
			int step = j;
			double current;
			double next;
			while(found==0){
				if(j==DistLength-1) {
					quantiles[i] = 1;
					found = 1;
				}
				else if(j==0){
					current = ((Double)modVals.get(j)).doubleValue();
					if(Scores[i] == current || Scores[i] > current){
						quantiles[i] = ((Double)modDist.get(j)).doubleValue();
						found = 1;
					}else{
						quantiles[i] = 0;
						found = 1;
					}
				}else{
					current = ((Double)modVals.get(j)).doubleValue();
					next = ((Double)modVals.get(j+1)).doubleValue();
					if(Scores[i] == current || (Scores[i] > current && Scores[i] < next)){
						quantiles[i] = ((Double)modDist.get(j)).doubleValue();
						found = 1;
					}else if(Scores[i] <current){
						step = Math.max(step/2,1);
						j = j - step;
					}else {
						step = Math.max(step/2,1);
						j = j + step;
					}
				}
			}
		}
		return quantiles;
	}
	
	public static String[][] getSignatures(String type){
		Vector numbers = new Vector();
		Vector forward = new Vector();
		Vector reversed = new Vector();
		
		char fsep = File.separatorChar;
		String curDir = System.getProperty("user.dir");
		String sigFN = curDir + fsep + "Models" + fsep + type + "_" + "Signatures" + ".txt";
		File distFile = new File(sigFN);
		
		try {
			FileReader inStream = new FileReader(distFile);
			BufferedReader in = new BufferedReader(inStream);
			String lineS;
			String numS;
			String sigsS;
			String forwS;
			String revS;
			int BreakP;
			while((lineS = in.readLine()) != null){
				BreakP = lineS.indexOf("\t");
				numS = lineS.substring(0, BreakP);
				sigsS = lineS.substring(BreakP+1);
				BreakP = sigsS.indexOf("\t");
				sigsS = sigsS.substring(BreakP+1);
				BreakP = sigsS.indexOf("\t");
				forwS = sigsS.substring(0, BreakP);
				revS = sigsS.substring(BreakP+1);
				numbers.add(numS);
				forward.add(forwS);
				reversed.add(revS);
			}
		} catch(Exception e) {
			System.out.println("Error reading file " + sigFN + "\n"+e);
		}
		int numMods = numbers.size();
		String[][] Signatures = new String[3][numMods];
		for(int i = 0; i < numMods; i++){
			Signatures[0][i] = (String)numbers.get(i).toString();
			Signatures[1][i] = (String)forward.get(i).toString();
			Signatures[2][i] = (String)reversed.get(i).toString();
		}
		return Signatures;
	}
	
	public static HashMap<String,MotifGroup> loadMotifGroups(String modelList, String modelType){
		char fsep = File.separatorChar;
		// String modelFolder = folder + fsep + modelType + "_models";
		// 2013-11-05 CLZ The user directly specifies the folder where the models will be found

		File f = new File(modelList);
		String modelFolder = f.getParent();
		
		HashMap<String,MotifGroup> Motifs = new HashMap<String,MotifGroup>();
		
		String listFile;
		listFile = modelList;
		String hashFileName = modelFolder + fsep + "all.grps";
		File hashFile = new File(hashFileName);
		if(hashFile.exists()){
			System.out.println("Reading existing java object file");
			try{
				FileInputStream fis = new FileInputStream(hashFileName);
				ObjectInputStream in = new ObjectInputStream(fis);
				Motifs = (HashMap<String,MotifGroup>)in.readObject();
				in.close();
			}catch (Exception e1){
				System.out.println("Error reading java object file " + hashFileName + "\n"+e1);
			}
			return Motifs;
		}else{
			System.out.println("Reading group files and writing java object file");
			try{
				FileReader inStream = new FileReader(listFile);
				BufferedReader in = new BufferedReader(inStream);
				String lineS;
				while((lineS = in.readLine()) != null){
					lineS = lineS.substring(0, lineS.length()-10);
					MotifGroup group = new MotifGroup(modelFolder,modelType,lineS);
					Motifs.put(lineS, group);
				}
				in.close();
			}
			catch (Exception e1){
				System.out.println("Error reading files"+"\n"+e1);
			}
			try{
				FileOutputStream fos = new FileOutputStream(hashFileName);
				ObjectOutputStream out = new ObjectOutputStream(fos);
				out.writeObject(Motifs);
				out.close();
			}catch (Exception e1){
				System.out.println("Error writing hashTable file"+"\n"+e1);
			}
			return Motifs;
		}
	}
	//getQuantiles for internal use by JAR3DMatlab or WebJAR3D functions,
	//overloaded for new file system
	public static double[] getQuantilesA(double[] Scores, MotifGroup group) {
		int n = Scores.length;
		double[] quantiles = new double[n];
		Vector modDist = new Vector();
		Vector modVals = new Vector();
		
		char fsep = File.separatorChar;

			String lineS;
			String ValueS;
			String DistS;
			double Value;
			double Dist;
			int BreakP;
			
			StringTokenizer st = new StringTokenizer(group.Distribution,"\n");

			lineS = st.nextToken();
			
			while(st.hasMoreTokens())
			{
				BreakP = lineS.indexOf(" ");
				ValueS = lineS.substring(0, BreakP);
				DistS = lineS.substring(BreakP+1);
				Dist = Double.parseDouble(DistS);
				modDist.add(Dist);
				Value = Double.parseDouble(ValueS);
				modVals.add(Value);
				lineS = st.nextToken();
			}
		
		int DistLength = modVals.size();
		
		for(int i = 0; i < n; i++){
			int found = 0;
			int j = DistLength/2;
			int step = j;
			double current;
			double next;
			while(found==0){
				if(j==DistLength-1) {
					quantiles[i] = 1;
					found = 1;
				}
				else if(j==0){
					current = ((Double)modVals.get(j)).doubleValue();
					if(Scores[i] == current || Scores[i] > current){
						quantiles[i] = ((Double)modDist.get(j)).doubleValue();
						found = 1;
					}else{
						quantiles[i] = 0;
						found = 1;
					}
				}else{
					current = ((Double)modVals.get(j)).doubleValue();
					next = ((Double)modVals.get(j+1)).doubleValue();
					if(Scores[i] == current || (Scores[i] > current && Scores[i] < next)){
						quantiles[i] = ((Double)modDist.get(j)).doubleValue();
						found = 1;
					}else if(Scores[i] <current){
						step = Math.max(step/2,1);
						j = j - step;
					}else {
						step = Math.max(step/2,1);
						j = j + step;
					}
				}
			}
		}
		return quantiles;
	}
	
	
}