package edu.bgsu.rna.jar3d.results;

import edu.bgsu.rna.jar3d.Sequence;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;

/**
 * This interface contains the functions for a SequenceResult. Sequence results
 * are the result of running one model over one sequence. This produces
 * interface works to link back to the Sequence run, and the Loop it was part
 * of.
 */
public interface SequenceResult {

    /**
     * The Query this Sequence was part of
     *
     * @return The Query.
     */
	public Query query();

    /**
     * The Sequence this sequence result came from.
     *
     * @return The Sequence.
     */
    public Sequence sequence();

    /**
     * The Loop this SequenceResult came from.
     *
     * @return The Loop.
     */
    public Loop loop();

    /**
     * The ID of the Loop this belongs to.
     *
     * @return The id.
     */
	public long loopId();

    /**
     * Get the query Id this Sequence result is part of.
     *
     * @return The query id.
     */
	public String queryId();
	
	/**
     * Get the correspondences for this sequences.
     *
     * @return The correspondences.
     */
	public String correspondences();

    /**
     * Get the sequence id of the sequence this belongs to.
     *
     * @return The sequence id.
     */
	public int sequenceId();

    /**
     * Get the motif Id of the motif that this was scored against.
     *
     * @return The motif id.
     */
	public String motifId();

    /**
     * Get the score of this result.
     *
     * @return The score.
     */
	public double score();

    /**
     * Get the edit interior distance of this result.
     *
     * @return The edit distance.
     */
	public int InteriorEditDistance();
	
	/**
     * Get the full edit interior distance of this result.
     *
     * @return The edit distance.
     */
	public int FullEditDistance();

    /**
     * Return 0 for HL, 0 or 1 for IL, 0, 1, or 2 for 3WJ, etc., depending on what rotation fits the model best
     *
     * @return true if rotated.
     */
	public int bestRotation();

    /**
     * Set the loop result this belongs to.
     *
     * @param result The LoopResult.
     */
	public void setLoopResult(LoopResult result);

    /**
     * Get the Loop result this belongs to.
     *
     * @return The loop result.
     */
	public LoopResult getLoopResult();
	
	/**
     * Get a boolean indicating if this passes the cutoff.
     *
     * @return The loop result.
     */
	public boolean cutoff();
	
	/**
     * Get a double indicating if the cutoff score.
     *
     * @return The loop result.
     */
	public double cutoffscore();
}
