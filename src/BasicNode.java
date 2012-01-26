
/**
 * This is a wrapper class for basic (non-branching) nodes
 * @author meg pirrung
 *
 */
class BasicNode extends Node {
   
	double maxLIns;
    double maxRIns;

	public BasicNode(Node prev, String type, int lI, int rI) {
		super(prev, type, lI, rI);
	}


	String getParams() {
		return null;
	}
	
}
