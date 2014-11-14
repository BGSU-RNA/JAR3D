package edu.bgsu.rna.jar3d.io.writers;

import java.util.List;

import edu.bgsu.rna.jar3d.results.LoopResult;
import edu.bgsu.rna.jar3d.results.SequenceResult;
import edu.bgsu.rna.jar3d.query.Query;

/**
 * This defines the interface for saving results.
 */
public interface ResultSaver {

	public void writeHeader() throws SaveFailed;

	public void save(List<LoopResult> results) throws SaveFailed;

	public void save(LoopResult results) throws SaveFailed;

	public void markFailure(String queryId) throws SaveFailed;

	public void cleanUp() throws SaveFailed;
	
	public void markAllDone(String queryId) throws SaveFailed;

	public void saveSequenceResults(List<SequenceResult> results, Query query) throws SaveFailed;
}
