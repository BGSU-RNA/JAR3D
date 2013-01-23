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

    /** Edit distance of this result. */
	private int editDistance;

    /** True if this was run in the rotated direction. */
	private boolean rotation;

    /** Loop result this belongs to. */
	private LoopResult result;

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
	public BasicSequenceResult(LoopResult result, Sequence sequence, double score, double percentile, int editDistance, boolean rotation) {
		this.result = result;
        this.sequence = sequence;
		this.score = score;
		this.percentile = percentile;
		this.editDistance = editDistance;
		this.rotation = rotation;
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
	public BasicSequenceResult(Sequence sequence, double score, double percentile, int editDistance, boolean rotation) {
		this(null, sequence, score, percentile, editDistance, rotation);
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
}
