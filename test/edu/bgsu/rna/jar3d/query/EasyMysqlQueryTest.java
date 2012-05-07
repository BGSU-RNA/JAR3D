package edu.bgsu.rna.jar3d.query;

import static org.junit.Assert.*;

import java.sql.SQLException;

import org.junit.Before;
import org.junit.Test;

import edu.bgsu.rna.jar3d.query.MysqlQuery;
import edu.bgsu.rna.jar3d.query.Query;

public class EasyMysqlQueryTest {

	private Query query;

	@Before
	public void setUp() throws SQLException {
		query = new MysqlQuery("joe", "hacker", "jdbc:mysql://localhost:3306/jar3d", "04279bc5-b6fe-4312-859e-35a68630d335");
	}

	@Test
	public void testEasyLoad() {
		assertEquals(query.getId(), "04279bc5-b6fe-4312-859e-35a68630d335");
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
		assertEquals(query.loopCount(), 1);
	}
	
	@Test
	public void testLoopSequence() {
		String sequence = query.getLoops().get(0).getSequences().get(0);
		assertEquals(sequence, "GAAAC*GGACC");
	}
}
