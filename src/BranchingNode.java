import java.util.*;

/**
 * This is a wrapper class for branching nodes
 * @author meg pirrung
 *
 */
public class BranchingNode extends Node {
	LinkedList children;
	int branches;
	
	/**
	 * This is the BranchingNode constructor
	 * @param prev contains a pointer to the previous node
	 * @param type contains the type of branching node
	 * @param b is the number of branches
	 * @param c is the linked list of child nodes (branches)
	 */
	public BranchingNode(Node prev, String type, int b, LinkedList c, int rI, int lI) {
		super(prev, type, rI, lI);
		children = c;
		branches = b;
	}

	
	String getParams() {
		return null;
	}
	
	public String generate()
  	{        
		return "";
  	}

}
