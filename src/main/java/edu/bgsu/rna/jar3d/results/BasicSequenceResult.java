package edu.bgsu.rna.jar3d.results;

import edu.bgsu.rna.jar3d.Sequence;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;

public class BasicSequenceResult implements SequenceResult {

	private double score;

	private double percentile;

	private int editDistance;

	private boolean rotation;

	private String sequenceId;

	private String motifId;

	private LoopResult result;

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
	public BasicSequenceResult(LoopResult result, double score, double percentile, int editDistance, boolean rotation) {
		this.result = result;
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
	public BasicSequenceResult(double score, double percentile, int editDistance, boolean rotation) {
		this(null, score, percentile, editDistance, rotation);
	}

    @Override
	public double score() {
		return score;
	}

    @Override
	public double percentile() {
		return percentile;
	}

    @Override
	public int editDistance() {
		return editDistance;
	}

    @Override
	public boolean isRotated() {
		return rotation;
	}

    @Override
	public Query query() {
		return result.getQuery();
	}

    @Override
	public long loopId() {
		return result.getLoop().getId();
	}

    @Override
	public String queryId() {
		return query().getId();
	}

    @Override
	public String sequenceId() {
		return sequenceId;
	}

    @Override
	public String motifId() {
		return motifId;
	}

	@Override
	public void setLoopResult(LoopResult result) {
		this.result = result;
	}

	@Override
	public LoopResult getLoopResult() {
		return result;
	}

    @Override
    public Loop loop() {
        return getLoopResult().getLoop();
    }

    @Override
    public Sequence sequence() {
        return null;
    }
}
