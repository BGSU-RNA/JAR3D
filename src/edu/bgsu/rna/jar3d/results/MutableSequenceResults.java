package edu.bgsu.rna.jar3d.results;


public class MutableSequenceResults implements SequenceResult {

	private String groupId;
	
	private String sequence;
	
	private double score;
	
	private double percentile;
	
	private int editDistance;
	
	private boolean rotation;
	
	/**
	 * @param groupId the groupId to set
	 */
	public void setGroupId(String groupId) {
		this.groupId = groupId;
	}

	/**
	 * @param sequence the sequence to set
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence;
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

	public MutableSequenceResults(String groupId, String sequence, double score, double percentile,
			int editDistance, boolean rotation) {
		
		this.groupId = groupId;
		this.sequence = sequence;
		this.score = score;
		this.percentile = percentile;
		this.editDistance = editDistance;
		this.rotation = rotation;
	}
	
	public String groupId() {
		return groupId;
	}

	public String sequence() {
		return sequence;
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
}