package edu.bgsu.rna.jar3d.io.writers;

import java.util.List;

import edu.bgsu.rna.jar3d.results.LoopResult;

public abstract class AbstractResultsSaver implements ResultSaver {
	
	/**
	 * Save results for parsing a single model against a several loops. This will
	 * save the aggregate and and individual information. 
	 * 
	 * @param results The results of parsing a whole loop.
	 * @throws SaveFailed if any problem occurs. 
	 */
	public void save(List<LoopResult> results) throws SaveFailed {
		for(LoopResult result: results) {
			save(result);
		}
	}
}
