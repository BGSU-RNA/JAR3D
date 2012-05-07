package edu.bgsu.rna.jar3d.query;

import java.util.List;

public interface Query extends Iterable<Loop> {

	public String getId();
	
	public String getILSetName();
	
	public String getHLSetName();
	
	public String modelType();
	
	public boolean onlyStructured();
	
	public List<Loop> getLoops();
	
	public Loop getLoop(int index);
	
	public int loopCount();
}
