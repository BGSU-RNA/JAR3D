package edu.bgsu.rna.jar3d.query;

import static org.junit.Assert.*;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import edu.bgsu.rna.jar3d.io.loaders.QueryLoader;
import edu.bgsu.rna.jar3d.io.loaders.QueryLoadingFailed;
import edu.bgsu.rna.jar3d.io.writers.DBLoader;
import edu.bgsu.rna.jar3d.query.Query;


public class ComplexDBLoaderTest {

	private QueryLoader loader;
	private Query query;

	@Before
	public void setUp() throws QueryLoadingFailed, SQLException {
		loader = new DBLoader("joe", "hacker", "jdbc:mysql://localhost:3306/jar3d");
		query = loader.load("6d4f3d19-3940-48d9-883e-2c8c8e945625");
	}

	@Test
	public void testEasyLoad() {
		assertEquals(query.getId(), "6d4f3d19-3940-48d9-883e-2c8c8e945625");
	}
	
	@Test
	public void testILSet() {
		assertEquals(query.getILSetName(), "IL0.6");
	}
	
	@Test
	public void testHLSet() {
		assertEquals(query.getHLSetName(), "HL0.2");
	}
	
	@Test
	public void testModelType() {
		assertEquals(query.modelType(), "default");
	}
	
	@Test
	public void testStructured() {
		assertFalse(query.onlyStructured());
	}
	
	@Test
	public void testLoopCount() {
		assertEquals(query.loopCount(), 2);
	}
	
	@Test
	public void testLoopType() {
		assertEquals(query.getLoop(1).getTypeString(), "HL");
	}

	@Test
	public void testLoopSequences() {
		List<String> sequences = new ArrayList<String>();
		sequences.add("CUUUG");
		sequences.add("CUUUG");
		sequences.add("CUUUG");
		List<String> found = new ArrayList<String>();
		for(String sequence: query.getLoop(1).getSequenceStrings()) {
			found.add(sequence);
		}
		assertEquals(found, sequences);
	}
}
