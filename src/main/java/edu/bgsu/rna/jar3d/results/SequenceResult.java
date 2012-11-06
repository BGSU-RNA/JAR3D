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

	public String sequenceId();

	public String motifId();

	public double score();

	public double percentile();

	public int editDistance();

	public boolean isRotated();

	public void setLoopResult(LoopResult result);

	public LoopResult getLoopResult();
}
