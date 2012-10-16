package edu.bgsu.rna.jar3d.results;

import edu.bgsu.rna.jar3d.query.Query;

public abstract class AbstractSequenceResult implements SequenceResult {

	public Query query() {
		// TODO Auto-generated method stub
		return null;
	}

	public int loopId() {
		return 0;
	}

	public String queryId() {
		// TODO Auto-generated method stub
		return query().getId();
	}

	public String sequenceId() {
		return null;
	}

	public String motifId() {
		// TODO Auto-generated method stub
		return null;
	}

	public String sequence() {
		// TODO Auto-generated method stub
		return null;
	}

	public double score() {
		// TODO Auto-generated method stub
		return 0;
	}

	public double percentile() {
		// TODO Auto-generated method stub
		return 0;
	}

	public int editDistance() {
		// TODO Auto-generated method stub
		return 0;
	}

	public boolean isRotated() {
		// TODO Auto-generated method stub
		return false;
	}

}