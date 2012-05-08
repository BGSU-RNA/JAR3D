package edu.bgsu.rna.jar3d.query;

import java.util.Iterator;
import java.util.List;

public class BasicLoop implements Loop {

	private List<String> loops;
	
	private String type;
	
	public BasicLoop(List<String> loops, String type) {
		this.loops = loops;
		this.type = type;
	}

	public List<String> getSequences() {
		return loops;
	}

	public String getType() {
		return type;
	}

	public Iterator<String> iterator() {
		return loops.iterator();
	}
}