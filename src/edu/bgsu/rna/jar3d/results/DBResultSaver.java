package edu.bgsu.rna.jar3d.results;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.List;

public class DBResultSaver implements ResultsSaver {

	private final Connection connection;
	
    private PreparedStatement insertLoopResult;
    
    private PreparedStatement insertSequenceResult;
	
	public DBResultSaver(String username, String password, String db) throws SQLException {
        connection = DriverManager.getConnection(db, username, password);
        String loopResultSQL = "insert into bygroup (id, meanscore, meanpercentile, meaneditdist, medianscore, medianpercentile, medianeditdist, signature, rotation, groupnum) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?);";
        String sequenceResultSQL = "insert into bysequence (id, seqnum, sequence, score, percentile, editdist, rotation, groupnum) values(?, ?, ?, ?, ?, ?, ?, ?);";
        insertLoopResult = connection.prepareStatement(loopResultSQL);
        insertSequenceResult = connection.prepareStatement(sequenceResultSQL);
	}
	
	public void save(LoopResult results, boolean status) {
		// TODO Do all sequences at once for speed up?
		// TODO When save a sequence update time.
		for(SequenceResult sequenceResult: results.sequenceResults()) {
			saveSequenceResult(sequenceResult, status);
		}

		try {
			insertLoopResult.setString(1, results.loopId());
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			int count = insertLoopResult.executeUpdate();
			
			if (count == 0) {
				// Throw something - no rows updated.
			}
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
//	"insert into bysequence (id, seqnum, sequence, score, percentile, editdist, rotation, groupnum) values('%s', %d, '%s', %f, %f, %d, %d, %s) ", id,m,currentL,groupScores[m-1],quants[m-1],minDist[m-1],reversed[index],groupName);

	private void saveSequenceResult(SequenceResult result, boolean status) {
		try {
			insertSequenceResult.setString(1, result.groupId());
			insertSequenceResult.setString(3, result.sequence());
			insertSequenceResult.setInt(6, result.editDistance());

		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			int count = insertSequenceResult.executeUpdate();
			
			if (count == 0) {
				// Throw something
			}
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void cleanUp() {
		try {
			connection.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	public void save(List<LoopResult> results) {
		for(LoopResult result: results) {
			save(result);
		}
	}

	public void save(LoopResult results) {
		save(results, true);
	}
}