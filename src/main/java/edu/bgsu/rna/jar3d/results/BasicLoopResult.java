package edu.bgsu.rna.jar3d.results;

import java.util.List;

import edu.bgsu.rna.jar3d.ArrayMath;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;

public final class BasicLoopResult implements LoopResult {

	private Loop loop;

	private final String modelId;

	private final String signature;

	private double meanInteriorEditDistance;

	private double medianInteriorEditDistance;
	
	private double meanScore;

	private double meanPercentile;

	private double meanFullEditDistance;
	
	private double medianFullEditDistance;

	private double medianPercentile;

	private double medianScore;

	private boolean rotation;

	private final List<SequenceResult> sequenceResults;

	private String correspondices;

    /**
     * Build a new BasicLoopResult. The given SequenceResults will have their
     * LoopResult set to this object.
     *
     * @param modelId The model id that was used.
     * @param rotation True if the loop was run in a rotated direction.
     * @param signature The base pairing signature of the loop used.
     * @param sequenceResults The individual SequenceResults.
     * @param correspondecies The correspondences between the model and
     * sequences.
     */
	public BasicLoopResult(String modelId, boolean rotation,
			String signature, List<SequenceResult> sequenceResults,
			String correspondecies) {
		this.sequenceResults = sequenceResults;
		this.modelId = modelId;
		this.rotation = rotation;
		this.signature = signature;
		this.correspondices = correspondecies;
		computeData();
	}

    /**
     * Compute the mean and median data. This goes through all Sequence results,
     * and computes the means and medians as well as setting the LoopResult of
     * each to this object.
     */
	private void computeData() {
		int numSeqs = sequenceResults.size();
		double[] scores = new double[numSeqs];
		double[] percentiles = new double[numSeqs];
		int[] fullEdDists = new int[numSeqs];
		int[] interiorEdDists = new int[numSeqs];

		
		for(int i = 0; i < numSeqs; i++){
			SequenceResult seqR = sequenceResults.get(i);
			seqR.setLoopResult(this);
			scores[i] = seqR.score();
			percentiles[i] = seqR.percentile();
			interiorEdDists[i] = seqR.InteriorEditDistance();
			fullEdDists[i] = seqR.FullEditDistance();
		}

		this.medianScore = ArrayMath.median(scores);
		this.meanScore = ArrayMath.mean(scores);
		this.meanPercentile = 100*ArrayMath.mean(percentiles);
		this.medianPercentile = 100*ArrayMath.median(percentiles);
		this.meanInteriorEditDistance = ArrayMath.mean(interiorEdDists);
		this.medianInteriorEditDistance = ArrayMath.median(interiorEdDists);
		this.meanFullEditDistance = ArrayMath.mean(fullEdDists);
		this.medianFullEditDistance = ArrayMath.median(fullEdDists);
	}

	public String modelId() {
		return modelId;
	}

	public double meanInteriorEditDistance() {
		return meanInteriorEditDistance;
	}
	
	public double meanFullEditDistance() {
		return meanFullEditDistance;
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

	public double medianInteriorEditDistance() {
		return medianInteriorEditDistance;
	}
	
	public double medianFullEditDistance() {
		return medianFullEditDistance;
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

	public Loop getLoop() {
		return loop;
	}

	public void setLoop(Loop loop) {
		this.loop = loop;
	}

	public Query getQuery() {
		return getLoop().getQuery();
	}

	public String getQueryID() {
		return getQuery().getId();
	}

	public List<SequenceResult> sequenceResults() {
		return sequenceResults;
	}
}
