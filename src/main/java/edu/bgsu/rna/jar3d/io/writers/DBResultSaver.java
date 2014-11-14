package edu.bgsu.rna.jar3d.io.writers;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.List;

import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.results.LoopResult;
import edu.bgsu.rna.jar3d.results.SequenceResult;

public class DBResultSaver extends AbstractResultsSaver {

	private final Connection connection;
	
    private PreparedStatement insertLoopResult;
    
    private PreparedStatement insertSequenceResult;
    
    private PreparedStatement updateSequenceQuery;
    
    private PreparedStatement updateLoopQuery;
    
    private PreparedStatement markLoopFailure;
    
    private PreparedStatement insertCorrespondenceResult;
    
    private PreparedStatement getSequenceResult;
    
    private Timestamp now;
	
	public DBResultSaver(String username, String password, String dbConnection) throws SQLException {
        connection = DriverManager.getConnection(dbConnection, username, password);
        String loopResultSQL = "insert into jar3d_results_by_loop (query_id, loop_id, motif_id, cutoff_percent, mean_cutoff_score, meanscore, meaninterioreditdist, meanfulleditdist, medianscore, medianinterioreditdist, medianfulleditdist, signature, rotation, correspondences) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);";
        String sequenceResultSQL = "insert into jar3d_results_by_loop_instance (query_id, seq_id, loop_id, cutoff, score, interioreditdist, fulleditdist, rotation, motif_id) values(?, ?, ?, ?, ?, ?, ?, ?, ?);";
        String correspondenceResultSQL = "insert into jar3d_node_position_results (result_id, seq_pos, node_id, node_pos, insert) values(? , ?, ?, ?, ?, ?)";
        String updateLoopSQL = "UPDATE jar3d_query_info SET status=1, time_completed=? WHERE query_id = ?;";
        String updateSequenceSQL = "UPDATE jar3d_query_sequences SET status=1, time_completed=? WHERE query_id = ? and seq_id = ? and loop_id = ?;";
        String failureSQL = "UPDATE jar3d_query_info SET status=-1, time_completed = ? WHERE query_id = ?;";
        String loopResultIDSQL = "SELECT * FROM jar3d_results_by_loop_instance WHERE query_id = ? AND seq_id = ? AND loop_id = ? AND motif_id = ?;";
        insertLoopResult = connection.prepareStatement(loopResultSQL);
        insertSequenceResult = connection.prepareStatement(sequenceResultSQL);
        insertCorrespondenceResult = connection.prepareStatement(correspondenceResultSQL);
        updateLoopQuery = connection.prepareStatement(updateLoopSQL);
        updateSequenceQuery = connection.prepareStatement(updateSequenceSQL);
        getSequenceResult = connection.prepareStatement(loopResultIDSQL);
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
	private void saveSequenceResult(SequenceResult result, Query query) throws SaveFailed {
		int cutoff=0;
		if(result.cutoff()){
			cutoff=1;
		}
		try {
			
			String seq_id = "0";
			if (!result.sequenceId().isEmpty()) {
				seq_id = result.sequenceId();
			}

			insertSequenceResult.setString(1, query.getId());
			insertSequenceResult.setString(2, seq_id);
			insertSequenceResult.setInt(3, (int)result.loopId());
			insertSequenceResult.setInt(4, cutoff);
			insertSequenceResult.setFloat(5, (float)Math.max(Math.min(result.score(), 9999),-9999));
			insertSequenceResult.setInt(6, Math.max(Math.min(result.InteriorEditDistance(),9999),-9999));
			insertSequenceResult.setInt(7, result.FullEditDistance());
			insertSequenceResult.setInt(8, result.bestRotation());
			insertSequenceResult.setString(9, result.motifId());
			updateSequenceQuery.setTimestamp(1, now);
			updateSequenceQuery.setString(2, query.getId());
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
	 * Save a the results for a set of sequences.
	 * 
	 * @param result The results of analyzing a single sequence.
	 * @throws SaveFailed If any problem occurs. s
	 */
	
	public void saveSequenceResults(List<SequenceResult> results, Query query) throws SaveFailed {
		for(SequenceResult sequenceResult: results) {
			saveSequenceResult(sequenceResult, query);
			saveSequenceCorrespondences(sequenceResult, query);
		}
	}
	
	/**
	 * Save a the correspondences for a sequences.
	 * 
	 * @param result The results of analyzing a single sequence.
	 * @throws SaveFailed If any problem occurs. s
	 */
	public void saveSequenceCorrespondences(SequenceResult result, Query query) throws SaveFailed {
		String query_id = query.getId();
		String seq_id = result.sequenceId();
		long loop_id = result.loopId();
		String motif_id = result.motifId();
		ResultSet seq_result;
		int res_id;
		try {
			getSequenceResult.setString(1, query_id);
			getSequenceResult.setString(2, seq_id);
			getSequenceResult.setLong(3, loop_id);
			getSequenceResult.setString(4, motif_id);
		} catch (SQLException e) {
			throw new SaveFailed("Could not generate sequence result info sql.", e);
		}
		try {
			seq_result = getSequenceResult.executeQuery();
			System.out.println(seq_result);
			res_id = seq_result.getInt("id");
		} catch (SQLException e) {
			throw new SaveFailed("Could not retrieve sequence result info.", e);
		}
		String[] corrs = result.correspondences().replaceAll("\\r", "").split("\\n");
		for(int i = 0; i < corrs.length; i++) {
			String line = corrs[i];
			if(line.contains("has")){continue;} // The "Has" lines are for information already in the DB
			String[] parts = line.split("_");
			int seq_pos = Integer.parseInt(parts[3]);
			int node = Integer.parseInt(parts[9]);
			int node_pos = Integer.parseInt(parts[11]);
			boolean insertion = false;
			if(line.contains("Insertion")){insertion = true;}
			try{
				insertCorrespondenceResult.setInt(1, res_id);
				insertCorrespondenceResult.setInt(2, seq_pos);
				insertCorrespondenceResult.setInt(3, node);
				insertCorrespondenceResult.setInt(4, node_pos);
				insertCorrespondenceResult.setBoolean(5, insertion);
			}catch (SQLException e) {
				throw new SaveFailed("Could not generate correpspondences SQL.", e);
			}
			try {
				int count = insertCorrespondenceResult.executeUpdate();
				insertCorrespondenceResult.executeUpdate();
				
				if (count == 0) {
					throw new SaveFailed("Saving should have updated at least 1 row.");
				}
			} catch (SQLException e) {
				throw new SaveFailed("Could not save correspondences.", e);
			}
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