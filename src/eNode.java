
/**
 * This node is used to signal the end of branches for alternative nodes
 * @author meg pirrung
 *
 */
public class eNode extends BasicNode {

	public eNode(Node child, int rI, int lI) {
		super(child, "eNode", rI, lI);
	}

	
	String getParams() {
		return null;
	}
}
