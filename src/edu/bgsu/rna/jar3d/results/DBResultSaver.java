package edu.bgsu.rna.jar3d.results;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.List;


public class DBResultSaver implements ResultsSaver {

	private final Connection connection;
	
    private PreparedStatement insertLoopResult;
    
    private PreparedStatement insertSequenceResult;
    
    private PreparedStatement updateSequenceQuery;
    
    private PreparedStatement updateLoopQuery;
    
    private PreparedStatement markLoopFailure;
    
    private Timestamp now;
	
	public DBResultSaver(String username, String password, String db) throws SQLException {
        connection = DriverManager.getConnection(db, username, password);
        String loopResultSQL = "insert into jar3d_results_by_loop (query_id, loop_id, motif_id, meanscore, meanpercentile, meaneditdist, medianscore, medianpercentile, medianeditdist, signature, rotation, correspondences) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);";
        String sequenceResultSQL = "insert into jar3d_results_by_loop_instance (query_id, seq_id, loop_id, score, percentile, editdist, rotation, motif_id) values(?, ?, ?, ?, ?, ?, ?, ?);";
        String updateLoopSQL = "UPDATE jar3d_query_info SET status=1, time_completed=? WHERE query_id = ?;";
        String updateSequenceSQL = "UPDATE jar3d_query_sequences SET status=1, time_completed=? WHERE query_id = ? and seq_id = ? and loop_id = ?;";
        String failureSQL = "UPDATE jar3d_query_info SET status=-1, time_completed = ? WHERE query_id = ?;";
        insertLoopResult = connection.prepareStatement(loopResultSQL);
        insertSequenceResult = connection.prepareStatement(sequenceResultSQL);
        updateLoopQuery = connection.prepareStatement(updateLoopSQL);
        updateSequenceQuery = connection.prepareStatement(updateSequenceSQL);
        markLoopFailure = connection.prepareStatement(failureSQL);
//        now = new Time(new Date().getTime());
      now = new Timestamp(System.currentTimeMillis());
	}
	
	/**
	 * Save the information for a single loop run against a single model. This
	 * save the aggregate and sequence level information. 
	 * 
	 * @param results The loop results.
	 */
	public void save(LoopResult results)  throws SaveFailed {
		// TODO Do all sequences at once for speed up?
		// TODO When save a sequence update time.
		for(SequenceResult sequenceResult: results.sequenceResults()) {
			saveSequenceResult(sequenceResult);
		}

		try {

			insertLoopResult.setString(1, results.queryId());
			insertLoopResult.setInt(2, results.loopId());
			insertLoopResult.setString(3, results.modelId());
			insertLoopResult.setFloat(4, (float)results.meanScore());
			insertLoopResult.setFloat(5, (float)results.meanPercentile());
			insertLoopResult.setFloat(6, (float)results.meanEditDistance());
			insertLoopResult.setFloat(7, (float)results.medianScore());
			insertLoopResult.setFloat(8, (float)results.medianPercentile());
			insertLoopResult.setFloat(9, (float)results.meanEditDistance());
			insertLoopResult.setString(10, results.signature());
			insertLoopResult.setInt(11, rotationInt(results.isRotated()));
			insertLoopResult.setString(12, "Intentially left empty.");
			updateLoopQuery.setTimestamp(1, now);
			updateLoopQuery.setString(2, results.queryId());
		} catch (SQLException e) {
			throw new SaveFailed("Could not generate loop sql.", e);
		}
		
		try {
			int count = insertLoopResult.executeUpdate();
			updateLoopQuery.executeUpdate();
			
			if (count == 0) {
				throw new SaveFailed("Update count wrong");
			}
		} catch (SQLException e) {
			throw new SaveFailed("Save failed.", e);
		}
	}

	/**
	 * Save a the results for a single sequence.
	 * 
	 * @param result The results of analyzing a single sequence.
	 * @throws SaveFailed If any problem occurs. s
	 */
	private void saveSequenceResult(SequenceResult result) throws SaveFailed {
		try {
			int rotated = 0;
			if (result.isRotated()) {
				rotated = 1;
			}

			insertSequenceResult.setString(1, result.queryId());
			insertSequenceResult.setString(2, result.sequenceId());
			insertSequenceResult.setInt(3, result.loopId());
			insertSequenceResult.setFloat(4, (float)result.score());
			insertSequenceResult.setFloat(5, (float)result.percentile());
			insertSequenceResult.setInt(6, result.editDistance());
			insertSequenceResult.setInt(7, rotated);
			insertSequenceResult.setString(8, result.motifId());
			updateSequenceQuery.setTimestamp(1, now);
			updateSequenceQuery.setString(2, result.queryId());
			updateSequenceQuery.setString(3, result.sequenceId());
			updateSequenceQuery.setInt(4, result.loopId());
		} catch (SQLException e) {
			throw new SaveFailed("Could not generate save sequence sql.", e);
		}
		
		try {
			int count = insertSequenceResult.executeUpdate();
			updateSequenceQuery.executeUpdate();
			
			if (count == 0) {
				throw new SaveFailed("Saving should have updated at least 1 row.");
			}
		} catch (SQLException e) {
			throw new SaveFailed("Could not save sequence.", e);
		}
	}
	
	/**
	 * Close the database connection.
	 */
	public void cleanUp() {
		try {
			connection.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}
	
	private int rotationInt(boolean rotation) {
		if (rotation) {
			return 1;
		}
		return 0;
	}

	/**
	 * Save results for parsing a single model against a several loop. This will
	 * save the aggregate and and individual information. 
	 * 
	 * @param results The results of parsing a whole loop.
	 * @throws SaveFailed if any problem occurs. 
	 */
	public void save(List<LoopResult> results) throws SaveFailed {
		for(LoopResult result: results) {
			save(result);
		}
	}

	/**
	 * Mark that we could not process the loop for some reason.
	 * 
	 * @param loop The loop to mark.
	 */
	public void markFailure(String queryId) throws SaveFailed {
		try {
			markLoopFailure.setTimestamp(1, now);
			markLoopFailure.setString(2, queryId);
			markLoopFailure.executeUpdate();
		} catch(Exception e) {
			throw new SaveFailed("Could not mark failure of: " + queryId, e);
		}
	}

	// No header to write.
	public void writeHeader() throws SaveFailed { }
}