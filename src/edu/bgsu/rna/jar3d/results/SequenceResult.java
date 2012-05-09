package edu.bgsu.rna.jar3d.results;

public interface SequenceResult {

	public String groupId();
	
	public String sequence();
	
	public double score();
	
	public double percentile();
	
	public int editDistance();
	
	public boolean isRotated();
}
