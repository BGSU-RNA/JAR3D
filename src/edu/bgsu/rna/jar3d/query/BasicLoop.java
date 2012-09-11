package edu.bgsu.rna.jar3d.query;

import java.util.Iterator;
import java.util.List;

public class BasicLoop implements Loop {

	private List<String> loops;
	
	private String type;
	
	private final long id;
	
	private final String name;
	
	public BasicLoop(String name, long id, List<String> loops, String type) {
		this.loops = loops;
		this.type = type;
		this.id = id;
		this.name = name;
	}
	
	public BasicLoop(long id, List<String> loops, String type) {
		this(null, id, loops, type);
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

	public long getId() {
		return id;
	}

	public String getName() {
		return name;
	}
}