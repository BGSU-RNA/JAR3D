package edu.bgsu.rna.jar3d.results;

import java.util.List;

import edu.bgsu.rna.jar3d.query.Query;

/**
 * This class represents the results of an entire query. This severs as a container of both the LoopResults and some
 * meta data.
 */
public interface QueryResult extends Iterable<LoopResult> {

	/**
	 * Get all LoopResults in this QueryResult.
	 * 
	 * @return The LoopResults.
	 */
	public List<LoopResult> getLoopResults();

	/**
	 * Set the query this Query result belongs to.
	 * 
	 * @param The query.
	 */
	public void setQuery(Query query);

	/**
	 * Get the query this QueryResult comes from.
	 * 
	 * @return The Query.
	 */
	public Query getQuery();

	/**
	 * Add a LoopResult to this QueryResult.
	 * 
	 * @param loopResult The result to add.
	 */
	public void addLoop(LoopResult loopResult);
}
