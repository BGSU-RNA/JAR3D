package edu.bgsu.rna.jar3d.results;

import java.util.List;

import edu.bgsu.rna.jar3d.ArrayMath;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;

public final class BasicLoopResult implements LoopResult {

	private Loop loop;

	private final int loopId;

	private final String modelId;

	private final String signature;

	public final String name;

	private double meanEditDistance;

	private double meanScore;

	private double meanPercentile;

	private double medianEditDistance;

	private double medianPercentile;

	private double medianScore;

	private boolean rotation;

	private final List<SequenceResult> sequenceResults;

	private String correspondices;

	public BasicLoopResult(int loopId, String modelId, boolean rotation, 
			String signature, List<SequenceResult> sequenceResults, 
			String correspondecies, String name) {
		this.sequenceResults = sequenceResults;
		this.modelId = modelId;
		this.loopId = loopId;
		this.rotation = rotation;
		this.signature = signature;
		this.correspondices = correspondecies;
		this.name = name;
		computeData();
	}
	
	private void computeData() {
		int numSeqs = sequenceResults.size();
		double[] scores = new double[numSeqs];
		double[] percentiles = new double[numSeqs];
		int[] edDists = new int[numSeqs];
		
		for(int i = 0; i < numSeqs; i++){
			SequenceResult seqR = sequenceResults.get(i);
			seqR.setLoopResult(this);
			scores[i] = seqR.score();
			percentiles[i] = seqR.percentile();
			edDists[i] = seqR.editDistance();
		}
		
		this.medianScore = ArrayMath.median(scores);
		this.meanScore = ArrayMath.mean(scores);
		this.meanPercentile = 100*ArrayMath.mean(percentiles);
		this.medianPercentile = 100*ArrayMath.median(percentiles);
		this.meanEditDistance = ArrayMath.mean(edDists);
		this.medianEditDistance = ArrayMath.median(edDists);
	}
	
	public int loopId() {
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

	public String correspondences() {
		return correspondices;
	}

	public String name() {
		return name;
	}

	@Override
	public Loop getLoop() {
		return loop;
	}

	@Override
	public void setLoop(Loop loop) {
		this.loop = loop;
	}

	@Override
	public Query getQuery() {
		return getLoop().getQuery();
	}

	@Override
	public List<SequenceResult> sequenceResults() {
		return sequenceResults;
	}
}
