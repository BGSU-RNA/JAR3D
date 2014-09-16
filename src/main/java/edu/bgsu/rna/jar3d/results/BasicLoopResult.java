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

	private double meanFullEditDistance;
	
	private double medianFullEditDistance;

	private double medianScore;
	
	private double meanCutoff;
	
	private double meanCutoffScore;

	private int rotation;

	private final List<SequenceResult> sequenceResults;

	private String correspondences;

    /**
     * Build a new BasicLoopResult. The given SequenceResults will have their
     * LoopResult set to this object.
     *
     * @param modelId The model id that was used.
     * @param rotation True if the loop was run in a rotated direction.
     * @param signature The base pairing signature of the loop used.
     * @param sequenceResults The individual SequenceResults.
     * @param correspondences The correspondences between the model and
     * sequences.
     */
	public BasicLoopResult(String modelId, int rotation,
			String signature, List<SequenceResult> sequenceResults,
			String correspondences) {
		this.sequenceResults = sequenceResults;
		this.modelId = modelId;
		this.rotation = rotation;
		this.signature = signature;
		this.correspondences = correspondences;
		computeData();
	}
	
	public BasicLoopResult(String modelId, int rotation,
			String signature, String correspondences, 
			double medianScore, double meanScore,
			double meanInteriorEditDistance, double medianInteriorEditDistance,
			double meanFullEditDistance, double medianFullEditDistance,
			double meanCutoff, double meanCutoffScore) {
		this.sequenceResults = null;
		this.modelId = modelId;
		this.rotation = rotation;
		this.signature = signature;
		this.correspondences = correspondences;
		this.medianScore = medianScore;
		this.meanScore = meanScore;
		this.meanInteriorEditDistance = meanInteriorEditDistance;
		this.medianInteriorEditDistance = medianInteriorEditDistance;
		this.meanFullEditDistance = meanFullEditDistance;
		this.medianFullEditDistance = medianFullEditDistance;
		this.meanCutoff = meanCutoff;
		this.meanCutoffScore = meanCutoffScore;
	}

    /**
     * Compute the mean and median data. This goes through all Sequence results
     * and computes the means and medians as well as setting the LoopResult of
     * each to this object.
     */
	private void computeData() {
		int numSeqs = sequenceResults.size();
		double[] scores = new double[numSeqs];
		double[] cutoffscores = new double[numSeqs];
		int[] fullEdDists = new int[numSeqs];
		int[] interiorEdDists = new int[numSeqs];
		boolean[] cutoffs = new boolean[numSeqs];

		
		for(int i = 0; i < numSeqs; i++){
			SequenceResult seqR = sequenceResults.get(i);
			seqR.setLoopResult(this);
			scores[i] = seqR.score();
			interiorEdDists[i] = seqR.InteriorEditDistance();
			fullEdDists[i] = seqR.FullEditDistance();
			cutoffs[i] = seqR.cutoff();
			cutoffscores[i] = seqR.cutoffscore();
		}

		this.medianScore = ArrayMath.median(scores);
		this.meanScore = ArrayMath.mean(scores);
		this.meanInteriorEditDistance = ArrayMath.mean(interiorEdDists);
		this.medianInteriorEditDistance = ArrayMath.median(interiorEdDists);
		this.meanFullEditDistance = ArrayMath.mean(fullEdDists);
		this.medianFullEditDistance = ArrayMath.median(fullEdDists);
		this.meanCutoff = ArrayMath.mean(cutoffs);
		this.meanCutoffScore = ArrayMath.mean(cutoffscores);
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

	public double medianScore() {
		return medianScore;
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

	public int bestRotation() {
		return rotation;
	}

	public String correspondences() {
		return correspondences;
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
	
	public double meanCutoff(){
		return meanCutoff;
	}
	
	public double meanCutoffScore(){
		return meanCutoffScore;
	}
}
