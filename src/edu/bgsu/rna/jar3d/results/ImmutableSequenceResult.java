package edu.bgsu.rna.jar3d.results;

public class ImmutableSequenceResult implements SequenceResult {

	private final String groupId;
	
	private final String sequence;
	
	private final double score;
	
	private final double percentile;
	
	private final int editDistance;
	
	private final boolean rotation;
	
	public ImmutableSequenceResult(String groupId, String sequence, double score, double percentile,
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
