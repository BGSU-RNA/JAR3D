package edu.bgsu.rna.jar3d.io.writers;

public class SaveFailed extends Exception {

	private static final long serialVersionUID = -4199978302576212482L;

	public SaveFailed() {
		super();
	}
	
	public SaveFailed(String msg) {
		super(msg);
	}
	
	public SaveFailed(String msg, Throwable t) {
		super(msg, t);
	}
	
	public SaveFailed(Throwable t) {
		super(t);
	}
}
