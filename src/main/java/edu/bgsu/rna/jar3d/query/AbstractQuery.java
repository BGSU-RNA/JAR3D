package edu.bgsu.rna.jar3d.query;

import java.util.Iterator;

import edu.bgsu.rna.jar3d.loop.Loop;

public abstract class AbstractQuery implements Query {

	public Iterator<Loop> iterator() {
		return getLoops().iterator();
	}

	public Loop getLoop(int index) {
		return getLoops().get(index);
	}

	public int loopCount() {
		return getLoops().size();
	}
}