package edu.bgsu.rna.jar3d.query;

import java.util.List;

import edu.bgsu.rna.jar3d.loop.Loop;

/**
 * This class represents a query that cannot be altered once built.
 */
public class ImmutableQuery extends AbstractQuery {

	/** The loops that in this query. */
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

	@Override
	public String modelType() {
		return modelType;
	}

	@Override
	public boolean onlyStructured() {
		return onlyStructured;
	}

	@Override
	public List<Loop> getLoops() {
		return loops;
	}

	@Override
	public String getId() {
		return id;
	}

	@Override
	public String getILSetName() {
		return ilSet;
	}

	@Override
	public String getHLSetName() {
		return hlSet;
	}
}