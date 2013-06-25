package edu.bgsu.rna.jar3d;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.ArrayList;
import java.util.Vector;

import edu.bgsu.rna.jar3d.SimpleAlign;
import edu.bgsu.rna.jar3d.Sequence;

public class SimpleAlignTest {
	private List<Sequence> firstIL;
	private List<Sequence> secondIL;
	private List<Sequence> firstHL;
	private List<Sequence> secondHL;
	
	private List<Sequence> firstILsmall;
	private List<Sequence> secondILsmall;
	private List<Sequence> firstHLsmall;
	private List<Sequence> secondHLsmall;
	
	private List<Sequence> testIL;
	
    @Before
    public void setUp() {
    	
    	firstIL = new Vector<Sequence>();
    	secondIL = new Vector<Sequence>();
    	firstHL = new Vector<Sequence>();
    	secondHL = new Vector<Sequence>();
    	
    	firstILsmall = new Vector<Sequence>();
    	secondILsmall = new Vector<Sequence>();
    	firstHLsmall = new Vector<Sequence>();
    	secondHLsmall = new Vector<Sequence>();
    	
    	firstIL.add(new Sequence("header", "*"));
        secondIL.add(new Sequence("header", "*"));
        firstHL.add(new Sequence("header", "*"));
        secondHL.add(new Sequence("header", "*"));
        
        firstILsmall.add(new Sequence("header", "*"));
        secondILsmall.add(new Sequence("header", "*"));
        firstHLsmall.add(new Sequence("header", "*"));
        secondHLsmall.add(new Sequence("header", "*"));
    	
        firstIL.add(new Sequence("firstIL", "gaaag*cuuuc"));
        
        secondIL.add(new Sequence("firstIL", "gaaag*cuuuc"));
        secondIL.add(new Sequence("secondIL", "guaag*cuuuc"));
        secondIL.add(new Sequence("thirdIL", "gauag*cuuuc"));
        secondIL.add(new Sequence("fouthIL", "gaaug*cuuuc"));
        secondIL.add(new Sequence("fifthIL", "uaaag*cuuuc"));
        secondIL.add(new Sequence("sixthIL", "gaaau*cuuuc"));
        secondIL.add(new Sequence("seventhIL", "gaaag*cauuc"));
        secondIL.add(new Sequence("eigthIL", "gaaag*cuauc"));
        secondIL.add(new Sequence("ninthIL", "gaaag*cuuac"));
        secondIL.add(new Sequence("tenthIL", "gaaag*auuuc"));
        secondIL.add(new Sequence("eleventhIL", "gaaag*cuuua"));
        secondIL.add(new Sequence("twelthIL", "gaag*cuuuc"));
        secondIL.add(new Sequence("thirteenthIL", "gaaag*cuuc"));
        secondIL.add(new Sequence("fourteenthIL", "gaaaag*cuuuc"));
        secondIL.add(new Sequence("fifteenthIL", "gaaag*cuuuuc"));

        firstHL.add(new Sequence("firstHL", "aaaaa"));
        
        secondHL.add(new Sequence("firstHL", "aaaaa"));
        secondHL.add(new Sequence("secondHL", "uaaaa"));
        secondHL.add(new Sequence("thirdHL", "auaaa"));
        secondHL.add(new Sequence("forthHL", "aauaa"));
        secondHL.add(new Sequence("fifthHL", "aaaua"));
        secondHL.add(new Sequence("sixthHL", "aaaau"));
        secondHL.add(new Sequence("seventhHL", "aaaa"));
        secondHL.add(new Sequence("eighthHL", "aaaaaa"));
        
        firstILsmall.add(new Sequence("firstIL", "aa*uu"));
        
        secondILsmall.add(new Sequence("firstIL", "aa*uu"));
        secondILsmall.add(new Sequence("secondIL", "a*uu"));
        secondILsmall.add(new Sequence("thirdIL", "aa*u"));
        secondILsmall.add(new Sequence("fouthIL", "ga*uu"));
        secondILsmall.add(new Sequence("fifthIL", "ag*uu"));
        secondILsmall.add(new Sequence("sixthIL", "aa*gu"));
        secondILsmall.add(new Sequence("seventhIL", "aa*ug"));
        secondILsmall.add(new Sequence("eigthIL", "aaa*uu"));
        secondILsmall.add(new Sequence("ninthIL", "aa*uuu"));
        
        firstHLsmall.add(new Sequence("firstIL", "aa"));
        
        secondHLsmall.add(new Sequence("firstIL", "aa"));
        secondHLsmall.add(new Sequence("secondIL", "a"));
        secondHLsmall.add(new Sequence("thirdIL", "aaa"));
        secondHLsmall.add(new Sequence("fouthIL", "ga"));
        secondHLsmall.add(new Sequence("fifthIL", "ag"));
        

    	testIL = new Vector<Sequence>();
    	testIL.add(new Sequence("header", "*"));
    	testIL.add(new Sequence("CCGACCG*CCAGUACG", "*"));
    	
    }
    
    @Test
    public void testILEditDistFull() {
        int[][] editDistA = SimpleAlign.calcILEditDistances(firstIL, secondIL, false, false, false);
        int[][] answer = new int[][]{{0,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};
        assertArrayEquals(editDistA,answer);
    }
    
    @Test
    public void testILEditDistInt() {
        int[][] editDistA = SimpleAlign.calcILEditDistances(firstIL, secondIL, false, false, true);
        int[][] answer = new int[][]{{0,1,1,1,0,0,1,1,1,0,0,1,1,1,1}};
        assertArrayEquals(editDistA,answer);
    }
    
    @Test
    public void testHLEditDistFull() {
        int[][] editDistA = SimpleAlign.calcHLEditDistances(firstHL, secondHL, false, false);
        int[][] answer = new int[][]{{0,1,1,1,1,1,1,1}};
        assertArrayEquals(editDistA,answer);
    }
    
    @Test
    public void testHLEditDistInt() {
        int[][] editDistA = SimpleAlign.calcHLEditDistances(firstHL, secondHL, false, true);
        int[][] answer = new int[][]{{0,0,1,1,1,0,1,1}};
        assertArrayEquals(editDistA,answer);
    }
    
    @Test
    public void testILsmallEditDistFull() {
        int[][] editDistA = SimpleAlign.calcILEditDistances(firstILsmall, secondILsmall, false, false, false);
        int[][] answer = new int[][]{{0,1,1,1,1,1,1,1,1}};
        assertArrayEquals(editDistA,answer);
    }
    
    @Test
    public void testILsmallEditDistInt() {
        int[][] editDistA = SimpleAlign.calcILEditDistances(firstILsmall, secondILsmall, false, false, true);
        int[][] answer = new int[][]{{0,0,0,0,0,0,0,1,1}};
        assertArrayEquals(editDistA,answer);
    }
    
    @Test
    public void testHLsmallEditDistFull() {
        int[][] editDistA = SimpleAlign.calcHLEditDistances(firstHLsmall, secondHLsmall, false, false);
        int[][] answer = new int[][]{{0,1,1,1,1}};
        assertArrayEquals(editDistA,answer);
    }
    
    @Test
    public void testHLsmallEditDistInt() {
        int[][] editDistA = SimpleAlign.calcHLEditDistances(firstHLsmall, secondHLsmall, false, true);
        int[][] answer = new int[][]{{0,0,1,0,0}};
        assertArrayEquals(editDistA,answer);
    }
    
    @Test
    public void testILspecifcFull() {
        int[][] editDistA = SimpleAlign.calcILEditDistances(testIL, testIL, false, false, false);
        int[][] answer = new int[][]{{0}};
        assertArrayEquals(editDistA,answer);
    }
    
    @Test
    public void testILspecifcInt() {
        int[][] editDistA = SimpleAlign.calcILEditDistances(testIL, testIL, false, false, true);
        int[][] answer = new int[][]{{0}};
        assertArrayEquals(editDistA,answer);
    }
}
