package edu.bgsu.rna.jar3d.results;


import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class DBResultSaverTest {
	
	private ResultsSaver saver;
	
	private Connection connection;
	
	private Statement statement;
	
//	public MutableSequenceResults(String motifId, double score, double percentile,
//			int editDistance, boolean rotation) {
	
	private LoopResult generateResult() {
		String motifId = "ABob";
//		Query query = new
		int loopId = 0;
		MutableSequenceResults r1 = new MutableSequenceResults(motifId, 0, 0, 0, false, 0);
		r1.setLoopId(loopId);
		MutableSequenceResults r2 = new MutableSequenceResults(motifId, 0, 0, 0, false, 0);
		r2.setLoopId(loopId);
		List<SequenceResult> seqs = new ArrayList<SequenceResult>();
		seqs.add(r1);
		seqs.add(r2);
		LoopResult result = new ImmutableLoopResult(loopId, motifId, false, "sign", seqs, "corr");
		return result;
	}
	
	@Before
	public void setUp() throws Exception {
		String db = "jdbc:mysql://localhost:3306/jar3d";
		String user = "joe";
		String pass = "hacker";
		saver = new DBResultSaver(user, pass, db);
		LoopResult result = generateResult();
		saver.save(result);
		connection = DriverManager.getConnection(db, user, pass);
        statement = connection.createStatement();
	}

	@After
	public void tearDown() throws Exception {
		saver.cleanUp();
		statement.executeUpdate("DELETE FROM jar3d_loop_results");
	}

	@Test
	public void testInsertLoopResults() {
		
	}

	@Test
	public void testInsertSequenceResults() {
		
	}

	@Test
	public void testUpdateQueryInfo() {
		
	}
	
	@Test
	public void testUpdateQuerySequences() {
		
	}
	
	@Test
	public void testUpdateFailedQueryInfo() {
		
	}
	
	@Test
	public void testUpdateFailedQuerySequences() {
		
	}
}
