package edu.bgsu.rna.jar3d.results;

import java.util.List;

public abstract class AbstractResultsSaver implements ResultsSaver {
	
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
