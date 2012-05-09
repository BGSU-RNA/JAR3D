package edu.bgsu.rna.jar3d.results;

import java.util.List;

public interface ResultsSaver {

	public void save(List<LoopResult> results);

	public void save(LoopResult results);
	
	public void cleanUp();
}
