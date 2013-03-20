package edu.bgsu.rna.jar3d;

import java.util.*; 

public class SimpleAlign {
	
	public static String[] getSeqStrings(List<Sequence> sData){
		//pulls sequence strings from a vector of Sequence objects
		
		int n = sData.size();
		String seqs[] = new String[n-1];
		for(int i = 1; i < n; i++){
			Sequence seq = (Sequence)sData.get(i);
			seqs[i-1] = (String)seq.letters;
		}
		return seqs;	
	}
	
	public static String[] SimpleAlignment(String seq1, String seq2){
		//does a very basic alignment between seq1 and seq2
		
		int n1 = seq1.length();
		int n2 = seq2.length();
		int[][] Scores = new int[n1+1][n2+1];
		int[][] TraceBack = new int[n1+1][n2+1];
		for(int i = 0; i <= n1; i ++){
			Scores[i][0] = 0;
			TraceBack[i][0] = 1;
		}
		for(int i = 0; i <= n2; i ++){
			Scores[0][i] = 0;
			TraceBack[0][i] = 2;
		}
		TraceBack[0][0] = -1;
		for(int i = 1; i <= n1; i ++){
			for(int j = 1; j <= n2; j ++){
				int i1 = Scores[i-1][j];
				int i2 = Scores[i][j-1];
				int M = 0;
				if(seq1.charAt(i)==seq2.charAt(j)) M = 1;
				int match = Scores[i-1][j-1] + M;
				int best = 1;
				if(i2>i1) best=2;
				if(match>best) best=3;
				if(best==1){
					TraceBack[i][j] = 1;
					Scores[i][j] = i1;
				}else if(best==2){
					TraceBack[i][j] = 2;
					Scores[i][j] = i2;
				}else {
					TraceBack[i][j] = 3;
					Scores[i][j] = match;
				}
			}
		}
		String[] aligned = new String[2];
		int i = n1;
		int j = n2;
		while(TraceBack[i][j] != -1){
			if(TraceBack[i][j]==1){
				aligned[0] = seq1.subSequence(i,i) + aligned[0];
				aligned[1] = "-" + aligned[1];
				i = i--;
			}else if(TraceBack[i][j]==2){
				aligned[1] = seq1.subSequence(i,i) + aligned[1];
				aligned[0] = "-" + aligned[0];
				j = j--;
			}else{
				aligned[0] = seq1.subSequence(i,i) + aligned[0];
				aligned[1] = seq1.subSequence(j,j) + aligned[1];
				i = i--;
				j = j--;
			}
		}
		return aligned;
	}
	
	public static int editDist(String seq1, String seq2){
		//calculates the Levenshtein distance between seq1 and seq2
		
		
		seq1 = seq1.toUpperCase().replaceAll("-", "");
		seq2 = seq2.toUpperCase().replaceAll("-", "");
		
		if(seq1.length()==0&seq2.length()==0) return 0;
		if(seq1.length()==0) 				  return seq2.length();
		if(seq2.length()==0) 				  return seq1.length();
		
		int n1 = seq1.length();
		int n2 = seq2.length();
		
		int[][] Scores = new int[n1+1][n2+1];
		for(int i = 0; i <= n1; i ++){
			Scores[i][0] = i;
		}
		for(int i = 1; i <= n2; i ++){
			Scores[0][i] = i;
		}
		for(int i = 1; i <= n1; i ++){
			for(int j = 1; j <= n2; j ++){
				if(seq1.charAt(i-1)==seq2.charAt(j-1))
					Scores[i][j] = Scores[i-1][j-1];
				else {
					int i1 = Scores[i-1][j];
					int i2 = Scores[i][j-1];
					int i3 = Scores[i-1][j-1];
					int min = i1;
					if(i2<i1) min = i2;
					if(i3<min) min = i3;
					Scores[i][j] = min+1;
				}
			}
		}
		return Scores[n1][n2];
	}
	
	public static int[][] calcILEditDistances(List<Sequence> sD1,List<Sequence> sD2,boolean reverse){
		return calcILEditDistances(sD1,sD2,reverse,false,true);
	}

	public static int[][] calcILEditDistances(List<Sequence> sD1,List<Sequence> sD2,boolean reverse, boolean Verbose){
		return calcILEditDistances(sD1,sD2,reverse,Verbose,true);
	}
	
	public static int[][] calcILEditDistances(List<Sequence> sD1,List<Sequence> sD2,boolean reverse,boolean Verbose,boolean Interior){
		//calculates the edit distance between every sequence in fasta file seqFile1
		//and every sequence in faste file seqFile2.  the first dim of the returned
		//2d array corresponds to the files in seqFile1, the second to seqFile2
		//this function assumes the fasta files involved are for internal loops
		//and use a '*' character to denote the break in strands
		
		//reverse sequences in seqFile1 if reverse is true
		if(reverse){
			sD1 = Alignment.reverse(sD1);
		}
		
		String[] seqs1 = getSeqStrings(sD1);
		String[] seqs2 = getSeqStrings(sD2);
		int n1 = seqs1.length;
		int n2 = seqs2.length;
		String[] seqs1left = new String[n1];
		String[] seqs1right = new String[n1];
		String[] seqs2left = new String[n2];
		String[] seqs2right = new String[n2];
		int[][] EdDists = new int[n1][n2];
		int breakpoint;

		//break up sequences to left and right substrings
		if(Interior == true){
			for(int i = 0; i < n1; i++){
				breakpoint = seqs1[i].indexOf("*");
				if(breakpoint > 2){
					seqs1left[i] = seqs1[i].substring(1, breakpoint-1);
				}
				else seqs1left[i] = "";
				if(breakpoint < seqs1[i].length()-3){
					seqs1right[i] = seqs1[i].substring(breakpoint+2,seqs1[i].length()-1);
				}
				else seqs1right[i] = "";
			}
			for(int i = 0; i < n2; i++){
				breakpoint = seqs2[i].indexOf("*");
				if(breakpoint > 2){
					seqs2left[i] = seqs2[i].substring(1, breakpoint-1);
				}
				else seqs2left[i] = "";
				if(breakpoint < seqs2[i].length()-3){
					seqs2right[i] = seqs2[i].substring(breakpoint+2,seqs2[i].length()-1);
				}
				else seqs2right[i] = "";
			}
		}else{
			for(int i = 0; i < n1; i++){
				breakpoint = seqs1[i].indexOf("*");
				seqs1left[i] = seqs1[i].substring(0, breakpoint);
				seqs1right[i] = seqs1[i].substring(breakpoint+1,seqs1[i].length());
			}
			for(int i = 0; i < n2; i++){
				breakpoint = seqs2[i].indexOf("*");
				seqs2left[i] = seqs2[i].substring(0, breakpoint);
				seqs2right[i] = seqs2[i].substring(breakpoint+1,seqs2[i].length());
			}			
		
		}
		int leftDist;
		int rightDist;
		for(int i = 0; i < n1; i ++){
			for(int j = 0; j < n2; j ++){
				leftDist = editDist(seqs1left[i],seqs2left[j]);
				rightDist = editDist(seqs1right[i],seqs2right[j]);
				EdDists[i][j]=leftDist+rightDist;
				if(Verbose){
					System.out.println(String.format("Finding minimum edit distance between %s and %s.", seqs1[i],seqs2[j]));
					System.out.println(String.format("Or, without flanking pairs,%s*%s and %s*%s.", seqs1left[i],seqs1right[i],seqs1left[j],seqs1right[j]));
					System.out.println(String.format("Left Edit Distance: %d",leftDist));
					System.out.println(String.format("Right Edit Distance: %d",rightDist));
					System.out.println(String.format("Total Edit Distance: %d",leftDist+rightDist));
				}
			}
		}
		return EdDists;
	}
	
	public static int[][] calcHLEditDistances(Vector<Sequence> sD1,Vector<Sequence> sD2){
		return calcHLEditDistances(sD1,sD2,false,true);
	}
	
	public static int[][] calcHLEditDistances(Vector<Sequence> sD1,Vector<Sequence> sD2, boolean Verbose){
		return calcHLEditDistances(sD1,sD2,Verbose,true);
	}
	
	public static int[][] calcHLEditDistances(Vector<Sequence> sD1,Vector<Sequence> sD2,boolean Verbose,boolean Interior){
		//calculates the edit distance between every sequence in fasta file seqFile1
		//and every sequence in faste file seqFile2.  the first dim of the returned
		//2d array corresponds to the files in seqFile1, the second to seqFile2
		//this function assumes the fasta files involved are for hairpin loops
	
		String[] seqs1 = getSeqStrings(sD1);
		String[] seqs2 = getSeqStrings(sD2);
		int n1 = seqs1.length;
		int n2 = seqs2.length;
		int[][] EdDists = new int[n1][n2];
		String seq1mf;
		String seq2mf;

		for(int i = 0; i < n1; i ++){
			for(int j = 0; j < n2; j ++){
				if(Interior == true){
					seq1mf = seqs1[i].substring(1, seqs1[i].length()-2);
					seq2mf = seqs2[j].substring(1, seqs2[j].length()-2);
				}else{
					seq1mf = seqs1[i];
					seq2mf = seqs2[j];
				}
				EdDists[i][j]=editDist(seq1mf,seq2mf);
				if(Verbose){
					System.out.println(String.format("Finding minimum edit distance between %s and %s.", seqs1[i],seqs2[j]));
					System.out.println(String.format("Or, without flanking pairs,%s and %s.", seq1mf,seq2mf));
					System.out.println(String.format("Total Edit Distance: %d",EdDists[i][j]));
				}
			}
		}
		return EdDists;
	}
}