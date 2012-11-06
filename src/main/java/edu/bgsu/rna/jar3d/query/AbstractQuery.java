package edu.bgsu.rna.jar3d.query;

import java.util.Iterator;

import edu.bgsu.rna.jar3d.loop.Loop;

public abstract class AbstractQuery implements Query {

	@Override
	public Iterator<Loop> iterator() {
		return getLoops().iterator();
	}

	@Override
	public Loop getLoop(int index) {
		return getLoops().get(index);
	}

	@Override
	public int loopCount() {
		return getLoops().size();
	}
}