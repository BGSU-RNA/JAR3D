package edu.bgsu.rna.jar3d.query;

import java.util.List;

public class ImmutableQuery extends AbstractQuery {

	private final List<Loop> loops;

	private final boolean onlyStructured;

	private final String ilSet;

	private final String hlSet;

	private final String modelType;

	private final String id;
	
	public ImmutableQuery(String id, List<Loop> loops, boolean onlyStructured, String ilSet, String hlSet, String modelType) {
		this.id = id;
		this.loops = loops;
		this.onlyStructured = onlyStructured;
		this.ilSet = ilSet;
		this.hlSet = hlSet;
		this.modelType = modelType;
	}
	
	public String modelType() {
		return modelType;
	}

	public boolean onlyStructured() {
		return onlyStructured;
	}

	public List<Loop> getLoops() {
		return loops;
	}

	public String getId() {
		return id;
	}

	public String getILSetName() {
		return ilSet;
	}

	public String getHLSetName() {
		return hlSet;
	}
}