package edu.bgsu.rna.jar3d.loop;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import edu.bgsu.rna.jar3d.Sequence;
import edu.bgsu.rna.jar3d.query.Query;

/**
 * BasicLoop is a simple implementation of the Loop interface.
 */
public class BasicLoop implements Loop {

	/** Sequences which are part of this loop. */
	private List<Sequence> sequences;

	/** The type of loop. */
	private final LoopType type;

	/** A numeric ID for this loop. */
	private final long id;

	/** A human readable name for this loop. */
	private final String name;

	/** The query this loop belongs to. */
	private Query query;

	/**
	 * Create a new BasicLoop.
	 *
	 * @param name Name of this loop.
	 * @param id Id of this loop.
	 * @param loops The Sequences which are part of this loop.
	 * @param type The type of loop.
	 */
	public BasicLoop(String name, long id, List<Sequence> loops, LoopType type) {
		this.sequences = loops;
		// Other parts assume there is a header sequence.
		sequences.add(0, new Sequence("", ""));
		this.type = type;
		this.id = id;
		this.name = name;
	}

	/**
	 * Create a new BasicLoop.
	 *
	 * @param name Name of this loop.
	 * @param id Id of this loop.
	 * @param loops The Sequences which are part of this loop.
	 * @param type The type of loop.
	 */
	public BasicLoop(String name, long id, List<Sequence> loops, String type) {
		this(name, id, loops, LoopType.fromString(type));
	}

	/**
	 * Create a new BasicLoop. The type of the loop is infered from the first sequence.
	 *
	 * @param name Name of this loop.
	 * @param id Id of this loop.
	 * @param loops The Sequences which are part of this loop.
	 */
	public BasicLoop(String name, long id, List<Sequence> loops) {
		this(name, id, loops, LoopType.fromSequence(loops.get(0)));
	}

	/**
	 * Create a new BasicLoop. The name is set to be the id as a String.
	 *
	 * @param id The id.
	 * @param loops The Sequences which are part of this Loop.
	 * @param type The type of Loop.
	 */
	public BasicLoop(long id, List<Sequence> loops, String type) {
		this(Long.valueOf(id).toString(), id, loops, type);
	}

	/**
	 * Create a new BasicLoop. The name is set to be the id as a String. The type is inferred from the first sequence.
	 *
	 * @param id The id.
	 * @param loops The Sequences which are part of this Loop.
	 */
	public BasicLoop(long id, List<Sequence> loops) {
		this(Long.valueOf(id).toString(), id, loops, LoopType.fromSequence(loops.get(0)));
	}

	public List<String> getSequenceStrings() {
		List<String> sequenceStrings = new ArrayList<String>();
		for(Sequence sequence: sequences) {
			sequenceStrings.add(sequence.getSequence().toUpperCase());
		}
		return sequenceStrings;
	}

	public String getTypeString() {
		return getLoopType().getShortName();
	}

	public Iterator<Sequence> iterator() {
		return getSequences().iterator();
	}

	public long getId() {
		return id;
	}

	public String getName() {
		return name;
	}

	public List<Sequence> getSequences() {
		return sequences;
	}

	public Query getQuery() {
		return query;
	}

	public void setQuery(Query query) {
		this.query = query;
	}

	public LoopType getLoopType() {
		return type;
	}

    public Loop reverse() {
        List<Sequence> reversed = new ArrayList<Sequence>();
        for (Sequence sequence : getSequences().subList(1, getSequences().size())) {
            reversed.add(sequence.reverse());
        }
        return new BasicLoop(getName(), getId(), reversed, getLoopType());
    }
}
