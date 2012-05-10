package edu.bgsu.rna.jar3d.results;

import edu.bgsu.rna.jar3d.query.Query;


public class MutableSequenceResults implements SequenceResult {
	
	private double score;
	
	private double percentile;
	
	private int editDistance;
	
	private boolean rotation;
	
	private Query query;
	
	private int loopId;
	
	private String sequenceId;
	
	private String motifId;

	/**
	 * Create a new MutableSequenceResults. This contains the information for
	 * running a single sequence against a single model.
	 * 
	 * @param motifId The ID of the motif used to build the model.
	 * @param score The score.
	 * @param percentile Percentile of the sequence.
	 * @param editDistance The edit distance.
	 * @param rotation True if the sequence was rotated relative to the model.
	 */
	public MutableSequenceResults(String motifId, double score, double percentile,
			int editDistance, boolean rotation) {
		
		this.motifId = motifId;
		this.score = score;
		this.percentile = percentile;
		this.editDistance = editDistance;
		this.rotation = rotation;
	}
	
	/**
	 * @param query the query to set
	 */
	public void setQuery(Query query) {
		this.query = query;
	}

	/**
	 * @param loopId the loopId to set
	 */
	public void setLoopId(int loopId) {
		this.loopId = loopId;
	}

	/**
	 * @param sequenceId the sequenceId to set
	 */
	public void setSequenceId(String sequenceId) {
		this.sequenceId = sequenceId;
	}

	/**
	 * @param motifId the motifId to set
	 */
	public void setMotifId(String motifId) {
		this.motifId = motifId;
	}
	
	/**
	 * @param score the score to set
	 */
	public void setScore(double score) {
		this.score = score;
	}

	/**
	 * @param percentile the percentile to set
	 */
	public void setPercentile(double percentile) {
		this.percentile = percentile;
	}

	/**
	 * @param editDistance the editDistance to set
	 */
	public void setEditDistance(int editDistance) {
		this.editDistance = editDistance;
	}

	/**
	 * @param rotation the rotation to set
	 */
	public void setRotation(boolean rotation) {
		this.rotation = rotation;
	}

	public double score() {
		return score;
	}

	public double percentile() {
		return percentile;
	}

	public int editDistance() {
		return editDistance;
	}

	public boolean isRotated() {
		return rotation;
	}

	public Query query() {
		return query;
	}

	public int loopId() {
		return loopId;
	}

	public String queryId() {
		return query().getId();
	}

	public String sequenceId() {
		return sequenceId;
	}

	public String motifId() {
		return motifId;
	}
}