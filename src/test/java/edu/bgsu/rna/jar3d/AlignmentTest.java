package edu.bgsu.rna.jar3d;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.ArrayList;

import edu.bgsu.rna.jar3d.Alignment;

public class AlignmentTest {

    private List<List<Double>> numbers;
    
    @Before
    public void setUp() {
        numbers = new ArrayList<List<Double>>();
        for (int i = 0; i < 3; i++) {
            List<Double> current = new ArrayList<Double>();
            for (int j = 0; j < 3; j++) {
                current.add(Double.valueOf(i + j));
            }
            numbers.add(current);
        }
    }

    @Test
    public void testVec2Array() {
        Double[][] val = Alignment.vec2array(numbers);
        Double[][] ans = {{0.0, 1.0, 2.0}, {1.0, 2.0, 3.0}, {2.0, 3.0, 4.0}};
        assertEquals(val, ans);
    }
}
