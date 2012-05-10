package edu.bgsu.rna.jar3d.results;

import java.util.List;

import edu.bgsu.rna.jar3d.query.Query;

public interface LoopResult {

	public Query query();
	
	public String queryId();
		
	public int loopId();

	public String modelId();

	public double meanEditDistance();

	public double meanScore();

	public double meanPercentile();

	public double medianScore();

	public double medianPercentile();

	public double medianEditDistance();

	public String signature();

	public boolean isRotated();
	
	public String correspondencies();

	public List<SequenceResult> sequenceResults();
}
