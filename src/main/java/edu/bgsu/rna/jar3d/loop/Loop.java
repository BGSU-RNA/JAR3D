package edu.bgsu.rna.jar3d.loop;

import java.util.List;

import edu.bgsu.rna.jar3d.Sequence;
import edu.bgsu.rna.jar3d.query.Query;

/**
 * A loop is a set of sequences which can be run against a model. This interface defines the basic functions of a loop.
 * It serves both as a container for the list of sequences, as well as a way to store some meta data about the loop.
 */
public interface Loop extends Iterable<Sequence> {

	/**
	 * Get a list of all sequences in this loop as strings.
	 * 
	 * @return The sequences as strings.
	 */
	public List<String> getSequenceStrings();

	/**
	 * Get the short name for the type of loop this is.
	 * 
	 * @return The type.
	 */
	public String getTypeString();

	/**
	 * Get the LoopType.
	 * 
	 * @return The loop type.
	 */
	public LoopType getLoopType();

	/**
	 * Get a numeric ID for this loop. 
	 * @return
	 */
	public long getId();

	/**
	 * Get a human readable name for this loop.
	 * 
	 * @return The name
	 */
	public String getName();

	/**
	 * Transform this loop into a list of Sequence objects.
	 * 
	 * @return This loop as a series of Sequences.
	 */
	public List<Sequence> getSequences();

	/**
	 * Get the query this loop belongs to.
	 * 
	 * @return The query.
	 */
	public Query getQuery();

	/**
	 * Set the Query this loop belongs to.
	 * 
	 * @param query The query.
	 */
	public void setQuery(Query query);

    /**
     * Create a new Loop with the strands in the loop reversed. The strands are
     * '*' seperated.
     *
     * @return The reversed loop.
     */
    public Loop reverse();
}
