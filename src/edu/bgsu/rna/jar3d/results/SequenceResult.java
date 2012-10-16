package edu.bgsu.rna.jar3d.results;

import edu.bgsu.rna.jar3d.query.Query;

public interface SequenceResult {

	public Query query();

	public int loopId();

	public String queryId();

	public String sequenceId();

	public String motifId();

	public double score();

	public double percentile();

	public int editDistance();

	public boolean isRotated();
}