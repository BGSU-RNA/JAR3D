package edu.bgsu.rna.jar3d.results;

import edu.bgsu.rna.jar3d.Sequence;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;

/**
 * BasicSequenceResult is the a simple implementation of the SequenceResult
 * interface.
 */
public class BasicSequenceResult implements SequenceResult {

    /** Score of this result. */
	private double score;

    /** Percentile of this result. */
	private double percentile;

    /** Interior edit distance of this result. */
	private int InteriorEditDistance;
	
    /** Full edit distance of this result. */
	private int FullEditDistance;

    /** 0 for HL, 0 or 1 for IL, 0, 1, or 2 for 3WJ, etc. */
	private int rotation;

    /** Loop result this belongs to. */
	private LoopResult result;

    /** Indicator if sequence meets cutoff requirements */
    private boolean cutoff;
	
    /** The sequence scored. */
    private final Sequence sequence;
    

	/**
	 * Create a new MutableSequenceResults. This contains the information for
	 * running a single sequence against a single model.
	 *
	 * @param result The LoopResult this sequence result belongs to.
	 * @param score The score.
	 * @param percentile Percentile of the sequence.
	 * @param editDistance The edit distance.
	 * @param rotation True if the sequence was rotated relative to the model.
	 */
	public BasicSequenceResult(LoopResult result, Sequence sequence, double score, double percentile, int InteriorEditDistance, int FullEditDistance, int rotation, boolean cutoff) {
		this.result = result;
        this.sequence = sequence;
		this.score = score;
		this.percentile = percentile;
		this.InteriorEditDistance = InteriorEditDistance;
		this.FullEditDistance = FullEditDistance;
		this.rotation = rotation;
		this.cutoff = cutoff;
	}

	/**
	 * Create a new MutableSequenceResults. This contains the information for
	 * running a single sequence against a single model.
	 *
	 * @param score The score.
	 * @param percentile Percentile of the sequence.
	 * @param editDistance The edit distance.
	 * @param rotation True if the sequence was rotated relative to the model.
	 */
	public BasicSequenceResult(Sequence sequence, double score, double percentile, int InteriorEditDistance, int FullEditDistance, int rotation, boolean cutoff) {
		this(null, sequence, score, percentile, InteriorEditDistance, FullEditDistance, rotation, cutoff);
	}

	public double score() {
		return score;
	}

	public double percentile() {
		return percentile;
	}

	public int FullEditDistance() {
		return FullEditDistance;
	}
	
	public int InteriorEditDistance() {
		return InteriorEditDistance;
	}

	public int bestRotation() {
		return rotation;
	}

	public Query query() {
        LoopResult result = getLoopResult();
        if (result != null) {
            return result.getQuery();
        }
		return null;
	}

	public long loopId() {
        Loop result = loop();
        if (result != null) {
            return result.getId();
        }
		return -1;
	}

	public String queryId() {
        Query query = query();
        if (query != null) {
            return query.getId();
        }
		return null;
	}

	public String sequenceId() {
        Sequence sequence = sequence();
        if (sequence != null) {
            return sequence.getId();
        }
		return null;
	}

	public String motifId() {
        LoopResult result = getLoopResult();
        if (result != null) {
            return result.modelId();
        }
		return null;
	}

	public void setLoopResult(LoopResult result) {
		this.result = result;
	}

	public LoopResult getLoopResult() {
		return result;
	}

    public Loop loop() {
        LoopResult result = getLoopResult();
        if (result != null) {
            return result.getLoop();
        }
		return null;
    }

    public Sequence sequence() {
        return sequence;
    }

	public boolean cutoff() {
		// TODO Auto-generated method stub
		return cutoff;
	}
}
