package edu.bgsu.rna.jar3d.query;

import java.util.List;

import edu.bgsu.rna.jar3d.loop.Loop;

/**
 * A query represents a set of loops to scan with models. This works as both a container for the loops and some meta
 * data about the query. For example, this contains information about what type of models should be used, the type
 * of internal and hairpin loop models that should be used and the model type.
 */
public interface Query extends Iterable<Loop> {

	/**
	 * The ID for this query.
	 * 
	 * @return The id
	 */
	public String getId();

	/**
	 * Get the name of the internal loop model set to use.
	 * 
	 * @return The set name.
	 */
	public String getILSetName();

	/**
	 * Get the name of the hairpin loop models to use.
	 * 
	 * @return The set name.
	 */
	public String getHLSetName();

	/**
	 * Get the type of model to use.
	 * 
	 * @return The model type.
	 */
	public String modelType();

	/**
	 * Return true if this query uses only the structured models.
	 * 
	 * @return If this query uses only structured models.
	 */
	public boolean onlyStructured();

	/**
	 * Get all loops in this Query.
	 * 
	 * @return The loops.
	 */
	public List<Loop> getLoops();

	/**
	 * Get the Loop at the given index in this Query.
	 * 
	 * @param index Index to get.
	 * @return The Loop.
	 */
	public Loop getLoop(int index);

	/**
	 * Get the number of loops in this query.
	 * 
	 * @return The count.
	 */
	public int loopCount();
}
