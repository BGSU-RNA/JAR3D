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
    
    /** Score based on cutoff values */
    private double cutoffscore;
    
    /** Correspondences string */
    private String correspondences;
	
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
	public BasicSequenceResult(LoopResult result, Sequence sequence, double score, int InteriorEditDistance, int FullEditDistance, int rotation, boolean cutoff, double cutoffscore, String correspondences) {
		this.result = result;
        this.sequence = sequence;
		this.score = score;
		this.InteriorEditDistance = InteriorEditDistance;
		this.FullEditDistance = FullEditDistance;
		this.rotation = rotation;
		this.cutoff = cutoff;
		this.cutoffscore = cutoffscore;
		this.correspondences = correspondences;
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
	public BasicSequenceResult(Sequence sequence, double score, int InteriorEditDistance, int FullEditDistance, int rotation, boolean cutoff, double cutoffscore) {
		this(null, sequence, score, InteriorEditDistance, FullEditDistance, rotation, cutoff,cutoffscore,"NA");
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
	public BasicSequenceResult(Sequence sequence, double score, int InteriorEditDistance, int FullEditDistance, int rotation, boolean cutoff, double cutoffscore, String correspondences) {
		this(null, sequence, score, InteriorEditDistance, FullEditDistance, rotation, cutoff,cutoffscore,correspondences);
	}

	public double score() {
		return score;
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

	public int sequenceId() {
        return sequence.getSeqId();
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
		return cutoff;
	}
	
	public double cutoffscore() {
		return cutoffscore;
	}

	public String correspondences() {
		return correspondences;
	}
}
