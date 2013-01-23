package edu.bgsu.rna.jar3d.loop;

import edu.bgsu.rna.jar3d.Sequence;

/**
 * This enum represents the possible types of loops this can process. Currently, only hairpin and internal are used.
 */
public enum LoopType {

	/** Type for all internal loops. */
	INTERNAL("internal", "IL", 2), 

	/** Type for all hairpin loops. */
	HAIRPIN("hairpin", "HL", 1), 

	/** Type for junction loops. */
	JUNCTION("junction", "JL", 3),

	/** Type for all other loop types. */
	UNKNOWN("unknown", "?L", -1);

	/** Short name for this loop type. */
	private final String shortName;

	/** Long name of this loop type. */
	private final String longName;

	/** The number of strands in this loop. */
	private int strands;

	/**
	 * Create a new LoopType.
	 * 
	 * @param longName The long name.
	 * @param shortName The short name.
	 * @param strands The number of strands.
	 */
	private LoopType(String longName, String shortName, int strands) {
		this.longName = longName;
		this.shortName = shortName;
		this.strands = strands;
	}

	/**
	 * Get the correct loop type from the given string. The string may be either a short or long name. If nothing 
	 * matches the given loop type then UNKNOWN is returned. 
	 * 
	 * @param type The type
	 * @return The matching LoopType.
	 */
	public static LoopType fromString(String type) {
		if (type.equalsIgnoreCase("IL") || type.equalsIgnoreCase("internal")) {
			return INTERNAL;
		}
		if (type.equalsIgnoreCase("HL") || type.equalsIgnoreCase("hairpin")) {
			return HAIRPIN;
		}
		if (type.equalsIgnoreCase("JL") || type.equalsIgnoreCase("junction")) {
			return JUNCTION;
		}
		return UNKNOWN;
	}

	/**
	 * Infer the loop type from a Sequence. If the sequence has no * then it is a hairpin. If it has one * then it is 
	 * an internal loop. Otherwise it is a junction loop.
	 * 
	 * @param sequence
	 * @return The loop type.
	 */
	public static LoopType fromSequence(Sequence sequence) {
		String seq = sequence.getSequence();
		int firstStar = seq.indexOf("*");
		int lastStar = seq.lastIndexOf("*");
		if (firstStar == -1) {
			return HAIRPIN;
		}
		if (firstStar == lastStar) {
			return INTERNAL;
		}
		return JUNCTION;
	}

	/**
	 * Get the long name of this loop.
	 * 
	 * @return The long name.
	 */
	public String getLongName() {
		return longName;
	}

	/**
	 * Get the short name of this loop.
	 * 
	 * @return The short name.
	 */
	public String getShortName() {
		return shortName;
	}

	/**
	 * Get the number of strands in this loop.
	 * 
	 * @return The number of strands.
	 */
	public int getStrandCount() {
		return strands;
	}

	/**
	 * Return a string representation of the loop. The short name is used as a representation.
	 * 
	 * @return A string of this loop.
	 */
	@Override
	public String toString() {
		return getShortName();
	}
}