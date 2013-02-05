package edu.bgsu.rna.jar3d.results;

import java.util.List;

import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.Query;

/**
 * This interface represents the result of parsing a loop against a single model. It links the result to the original
 * loop for ease of tracking.
 */
public interface LoopResult {

	/**
	 * Query this result belongs to.
	 * 
	 * @return The query
	 */
	public Query getQuery();

	/**
	 * Get the loop that was scanned.
	 * 
	 * @return The loop
	 */
	public Loop getLoop();

	/**
	 * Set the loop this result belongs to.
	 * 
	 * @param loop The loop.
	 */
	public void setLoop(Loop loop);

	/**
	 * ID of the model scanned.
	 * 
	 * @return The model id.
	 */
	public String modelId();

	/**
	 * Mean interior edit distance between sequences in the loop and the sequences in the model.
	 * 
	 * @return The edit distance.
	 */
	public double meanInteriorEditDistance();
	
	/**
	 * Mean full interior edit distance between sequences in the loop and the sequences in the model.
	 * 
	 * @return The edit distance.
	 */
	public double meanFullEditDistance();

	/**
	 * The mean score of sequences in this loop.
	 * 
	 * @return The score.
	 */
	public double meanScore();

	/**
	 * The mean percentile score of sequences in this loop.
	 * 
	 * @return The percentile.
	 */
	public double meanPercentile();

	/**
	 * The median score of sequences in this loop.
	 * 
	 * @return The score.
	 */
	public double medianScore();

	/**
	 * The median percentile score of sequences in this loop.
	 * 
	 * @return The percentile.
	 */
	public double medianPercentile();

	/**
	 * Median interior edit distance between sequences in the loop and the sequences in the model.
	 * 
	 * @return The edit distance.
	 */
	public double medianInteriorEditDistance();
	
	/**
	 * Median full edit distance between sequences in the loop and the sequences in the model.
	 * 
	 * @return The edit distance.
	 */
	public double medianFullEditDistance();

	/**
	 * This is the base pairing signature for the scanned model.
	 * 
	 * @return The signature
	 */
	public String signature();

	/**
	 * True if the scanned loop was rotated. As in was this loop run in the reverse direction against the model.
	 * 
	 * @return True if rotated.
	 */
	public boolean isRotated();

	/**
	 * Get the correspondences between the model and positions in this loop.
	 * 
	 * @return The correspondences.
	 */
	public String correspondences();

	/**
	 * A list of per sequence results. These are the results for each sequence in loop against the model used here.
	 * 
	 * @return A list of sequence level results.
	 */
	public List<SequenceResult> sequenceResults();
}
