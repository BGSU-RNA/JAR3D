package edu.bgsu.rna.jar3d.results;

import java.util.List;

public final class ImmutableLoopResult implements LoopResult {

	private final String loopId;
	
	private final String modelId;
	
	private final String signature;
	
	private double meanEditDistance;
	
	private  double meanScore;
	
	private double meanPercentile;
	
	private double medianEditDistance;
	
	private double medianPercentile;
	
	private double medianScore;
	
	private boolean rotation;
	
	private final List<SequenceResult> sequenceResults;
	
	public ImmutableLoopResult(String loopId, String modelId, boolean rotation, String signature, List<SequenceResult> sequenceResults) {
		this.sequenceResults = sequenceResults;
		this.modelId = modelId;
		this.loopId = loopId;
		this.rotation = rotation;
		this.signature = signature;
		computeData();
	}
	
	private void computeData() {
		
	}
	
	public String loopId() {
		return loopId;
	}

	public String modelId() {
		return modelId;
	}

	public double meanEditDistance() {
		return meanEditDistance;
	}

	public double meanScore() {
		return meanScore;
	}

	public double meanPercentile() {
		return meanPercentile;
	}

	public double medianScore() {
		return medianScore;
	}

	public double medianPercentile() {
		return medianPercentile;
	}

	public double medianEditDistance() {
		return medianEditDistance;
	}

	public String signature() {
		return signature;
	}

	public boolean isRotated() {
		return rotation;
	}

	public List<SequenceResult> sequenceResults() {
		return sequenceResults;
	}

}
