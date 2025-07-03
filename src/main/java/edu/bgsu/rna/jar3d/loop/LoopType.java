package edu.bgsu.rna.jar3d.loop;

import edu.bgsu.rna.jar3d.Sequence;

/**
 * This enum represents the possible types of loops this code can process.
 */
public enum LoopType {

	/** Type for all hairpin loops. */
	HAIRPIN("hairpin", "HL", 1),

	/** Type for all internal loops. */
	INTERNAL("internal", "IL", 2),

	/** Type for junction loops. */
	JUNCTION("junction", "JL", 0),

	/** Type for 3-way junction loops. */
	J3("junction", "J3", 3),

	/** Type for 4-way junction loops. */
	J4("junction", "J4", 4),

	/** Type for 4-way junction loops. */
	J5("junction", "J5", 5),

	/** Type for 4-way junction loops. */
	J6("junction", "J6", 6),

	/** Type for 4-way junction loops. */
	J7("junction", "J7", 7),

	/** Type for 4-way junction loops. */
	J8("junction", "J8", 8),

	/** Type for 4-way junction loops. */
	J9("junction", "J9", 9),

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
		System.out.println(type);
		if (type.equalsIgnoreCase("IL") || type.equalsIgnoreCase("internal")) {
			return INTERNAL;
		}
		if (type.equalsIgnoreCase("HL") || type.equalsIgnoreCase("hairpin")) {
			return HAIRPIN;
		}
		if (type.equalsIgnoreCase("J3")) {
			return J3;
		}
		if (type.equalsIgnoreCase("J4")) {
			return J4;
		}
		if (type.equalsIgnoreCase("J5")) {
			return J5;
		}
		if (type.equalsIgnoreCase("J6")) {
			return J6;
		}
		if (type.equalsIgnoreCase("J7")) {
			return J7;
		}
		if (type.equalsIgnoreCase("J8")) {
			return J8;
		}
		if (type.equalsIgnoreCase("J9")) {
			return J9;
		}
		if (type.equalsIgnoreCase("JL") || type.equalsIgnoreCase("junction")) {
			return JUNCTION;
		}
		return UNKNOWN;
	}

	/**
	 * Infer the loop type from a Sequence.
	 * Use the number of * characters to count the strands.
	 *
	 * @param sequence
	 * @return The loop type.
	 */
	public static LoopType fromSequence(Sequence sequence) {
		String seq = sequence.getSequence();

		String[] parts = seq.split("\\*");  // escape * with \\, since it's a regex
        int strand_count = parts.length;

		if (strand_count == 1) {
			return HAIRPIN;
		} else if (strand_count == 2) {
			return INTERNAL;
		} else if (strand_count == 3) {
			return J3;
		} else if (strand_count == 4) {
			return J4;
		} else if (strand_count == 5) {
			return J5;
		} else if (strand_count == 6) {
			return J6;
		} else if (strand_count == 7) {
			return J7;
		} else if (strand_count == 8) {
			return J8;
		} else if (strand_count == 9) {
			return J9;
		}

		return UNKNOWN;

		// int firstStar = seq.indexOf("*");
		// int lastStar = seq.lastIndexOf("*");
		// if (firstStar == -1) {
		// 	return HAIRPIN;
		// }
		// if (firstStar == lastStar) {
		// 	return INTERNAL;
		// }
		// return JUNCTION;
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