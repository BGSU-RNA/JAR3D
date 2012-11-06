package edu.bgsu.rna.jar3d;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.ArrayList;

import edu.bgsu.rna.jar3d.Sequence;

public class SequenceTest {

    private Sequence sequence;

    private Sequence hairpin;

    @Before
    public void setUp() {
        sequence = new Sequence("first", "aaa*ccc");
        hairpin = new Sequence("hair", "aaa");
    }

    @Test
    public void testInternalReversalSequence() {
        Sequence reversed = sequence.reverse();
        assertEquals(reversed.getSequence(), "ccc*aaa");
    }

    @Test
    public void testHairpinReveralSequence() {
        Sequence reversed  = hairpin.reverse();
        assertEquals(reversed.getSequence(), "aaa");
    }
}
