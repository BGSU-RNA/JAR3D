package edu.bgsu.rna.jar3d.loop;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.ArrayList;

import edu.bgsu.rna.jar3d.Sequence;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.loop.BasicLoop;
import edu.bgsu.rna.jar3d.loop.LoopType;

public class BasicLoopTest {

    private Loop loop;

    @Before
    public void setUp() {
        List<Sequence> sequences = new ArrayList<Sequence>();
        sequences.add(new Sequence("first", "aaa*ccc"));
        loop = new BasicLoop("a", 1L, sequences, LoopType.INTERNAL);
    }

    @Test
    public void testLoopReversalSequences() {
        Loop reversed = loop.reverse();
        List<Sequence> sequences = reversed.getSequences();
        assertEquals(sequences.get(1).getSequence(), "ccc*aaa");
    }

    @Test
    public void testLoopReversalName() {
        Loop reversed = loop.reverse();
        assertEquals(reversed.getName(), "a");
    }

    @Test
    public void testLoopReversalId() {
        Loop reversed = loop.reverse();
        assertEquals(reversed.getId(), 1L);
    }

    @Test
    public void testLoopReversalType() {
        Loop reversed = loop.reverse();
        assertEquals(reversed.getLoopType(), LoopType.INTERNAL);
    }
}
