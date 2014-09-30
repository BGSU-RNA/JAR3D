package edu.bgsu.rna.jar3d.io.writers;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Timestamp;

import edu.bgsu.rna.jar3d.results.LoopResult;
import edu.bgsu.rna.jar3d.results.SequenceResult;

public class DBResultSaver extends AbstractResultsSaver {

	private final Connection connection;
	
    private PreparedStatement insertLoopResult;
    
    private PreparedStatement insertSequenceResult;
    
    private PreparedStatement updateSequenceQuery;
    
    private PreparedStatement updateLoopQuery;
    
    private PreparedStatement markLoopFailure;
    
    private Timestamp now;
	
	public DBResultSaver(String username, String password, String dbConnection) throws SQLException {
        connection = DriverManager.getConnection(dbConnection, username, password);
        String loopResultSQL = "insert into jar3d_results_by_loop (query_id, loop_id, motif_id, cutoff_percent, mean_cutoff_score, meanscore, meaninterioreditdist, meanfulleditdist, medianscore, medianinterioreditdist, medianfulleditdist, signature, rotation, correspondences) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);";
        String sequenceResultSQL = "insert into jar3d_results_by_loop_instance (query_id, seq_id, loop_id, cutoff, score, interioreditdist, fulleditdist, rotation, motif_id) values(?, ?, ?, ?, ?, ?, ?, ?, ?);";
        String updateLoopSQL = "UPDATE jar3d_query_info SET status=1, time_completed=? WHERE query_id = ?;";
        String updateSequenceSQL = "UPDATE jar3d_query_sequences SET status=1, time_completed=? WHERE query_id = ? and seq_id = ? and loop_id = ?;";
        String failureSQL = "UPDATE jar3d_query_info SET status=-1, time_completed = ? WHERE query_id = ?;";
        insertLoopResult = connection.prepareStatement(loopResultSQL);
        insertSequenceResult = connection.prepareStatement(sequenceResultSQL);
        updateLoopQuery = connection.prepareStatement(updateLoopSQL);
        updateSequenceQuery = connection.prepareStatement(updateSequenceSQL);
        markLoopFailure = connection.prepareStatement(failureSQL);
        now = new Timestamp(System.currentTimeMillis());
	}
	
	public void markAllDone(String queryId) throws SaveFailed {
		try {
			updateLoopQuery.setTimestamp(1, now);
			updateLoopQuery.setString(2, queryId);
			updateLoopQuery.executeUpdate();
		} catch (SQLException e) {
			throw new SaveFailed("Could not mark query as done", e);
		}
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
		// Disabling sequence saving to help with speed.  When sequence pages are displayed,
		// we need to re-do the calculations to get the correspondences, so we can save sequence level
		// results then
/*		for(SequenceResult sequenceResult: results.sequenceResults()) {
			saveSequenceResult(sequenceResult);
		}
*/
		try {
			insertLoopResult.setString(1, results.getLoop().getQuery().getId());
			insertLoopResult.setInt(2, (int)results.getLoop().getId());
			insertLoopResult.setString(3, results.modelId());
			insertLoopResult.setFloat(4, (float)Math.max(Math.min(results.meanCutoff(),9999),-9999));
			insertLoopResult.setFloat(5, (float)Math.max(Math.min(results.meanCutoffScore(),9999),-9999));
			insertLoopResult.setFloat(6, (float)Math.max(Math.min(results.meanScore(),9999),-9999));
			insertLoopResult.setFloat(7, (float)Math.max(Math.min(results.meanInteriorEditDistance(),9999),-9999));
			insertLoopResult.setFloat(8, (float)Math.max(Math.min(results.meanFullEditDistance(),9999),-9999));
			insertLoopResult.setFloat(9, (float)Math.max(Math.min(results.medianScore(),9999),-9999));
			insertLoopResult.setFloat(10, (float)Math.max(Math.min(results.medianInteriorEditDistance(),9999),-9999));
			insertLoopResult.setFloat(11, (float)Math.max(Math.min(results.medianFullEditDistance(),9999),-9999));
			insertLoopResult.setString(12, results.signature());
			insertLoopResult.setInt(13, results.bestRotation());
			insertLoopResult.setString(14, "Intentially left empty.");
		} catch (SQLException e) {
			throw new SaveFailed("Could not generate loop sql.", e);
		}
		
		try {
			int count = insertLoopResult.executeUpdate();

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
		int cutoff=0;
		if(result.cutoff()){
			cutoff=1;
		}
		try {
			
			String seq_id = "0";
			if (!result.sequenceId().isEmpty()) {
				seq_id = result.sequenceId();
			}

			insertSequenceResult.setString(1, result.queryId());
			insertSequenceResult.setString(2, seq_id);
			insertSequenceResult.setInt(3, (int)result.loopId());
			insertSequenceResult.setInt(4, cutoff);
			insertSequenceResult.setFloat(5, (float)Math.max(Math.min(result.score(), 9999),-9999));
			insertSequenceResult.setInt(6, Math.max(Math.min(result.InteriorEditDistance(),9999),-9999));
			insertSequenceResult.setInt(7, result.FullEditDistance());
			insertSequenceResult.setInt(8, result.bestRotation());
			insertSequenceResult.setString(9, result.motifId());
			updateSequenceQuery.setTimestamp(1, now);
			updateSequenceQuery.setString(2, result.queryId());
			updateSequenceQuery.setString(3, result.sequenceId());
			updateSequenceQuery.setInt(4, (int)result.loopId());
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