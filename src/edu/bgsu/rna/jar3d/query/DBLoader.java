package edu.bgsu.rna.jar3d.query;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import java.util.ArrayList;
import java.util.List;

/**
 * This is a class to load query data from the database.
 */
public class DBLoader implements QueryLoader {
    private Connection connection;

    private PreparedStatement sqlForQueryInfo;

    private PreparedStatement sqlForLoops;
    
    private PreparedStatement updateQueryInfo;

    public DBLoader(String username, String password, String dbConnection) throws SQLException {
        connection = DriverManager.getConnection(dbConnection, username, password);
        String querySql = "SELECT group_set, model_type, structured_models_only FROM `jar3d_query_info` WHERE query_id = ?;";
        String loopSql = "SELECT loop_id, loop_sequence, loop_type FROM `jar3d_query_sequences` WHERE query_id = ? and loop_id = ? and status = 0;";
        String updateInfo = "UPDATE jar3d_query_info SET status=2 WHERE query_id = ?;";
        sqlForQueryInfo = connection.prepareStatement(querySql);
        sqlForLoops = connection.prepareStatement(loopSql);
        updateQueryInfo = connection.prepareStatement(updateInfo);
    }

    private Loop loadLoop(String queryId, int index) throws SQLException {
    	sqlForLoops.setString(1, queryId);
    	sqlForLoops.setLong(2, index);
    	ResultSet results = sqlForLoops.executeQuery();
    	long id = -1;
    	
    	String type = null;
    	List<String> sequences = new ArrayList<String>();
    	while (results.next()) {
    		String sequence = results.getString("loop_sequence");
    		type = results.getString("loop_type");
    		id = results.getLong("loop_id");
    		sequences.add(sequence);
    	}
    	results.close();

    	return new BasicLoop(id, sequences, type);
    }

    private List<Loop> loadLoops(String queryId) throws SQLException, QueryLoadingFailed {
        String loopCountSql = "SELECT MAX(loop_id) AS max FROM `jar3d_query_sequences` where query_id = ? and status = 0;";
        PreparedStatement sqlForLoopCount = connection.prepareStatement(loopCountSql);
        sqlForLoopCount.setString(1, queryId);
        
        List<Loop> loops = new ArrayList<Loop>();
        ResultSet result = sqlForLoopCount.executeQuery();
        boolean found = result.first();
        if (!found) {
        	throw new QueryLoadingFailed("Could not find any loops for query: " + queryId);
        }
        
        int loopCount = result.getInt("max") + 1;
        result.close();
        
        for (int i = 0; i < loopCount; i++) {
        	loops.add(loadLoop(queryId, i));
        }
        
        return loops;
    }

	public Query load(String queryId) throws QueryLoadingFailed {

        String modelType;
        String groupSet;
        String ilSet;
        String hlSet;
        List<Loop> loops;
        boolean onlyStructured;
        
		try {
	        sqlForQueryInfo.setString(1, queryId);
	        updateQueryInfo.setString(1, queryId);
	        
	        updateQueryInfo.executeUpdate();
	        ResultSet result = sqlForQueryInfo.executeQuery();
	        boolean found = result.first();
	        
	        if (!found) {
	        	throw new QueryLoadingFailed("Could not find query with id: " + queryId);
	        }
	        
	        modelType = result.getString("model_type");
	        groupSet = result.getString("group_set");
	        String[] parts = groupSet.split("/");
	        ilSet = parts[0];
	        hlSet = parts[1];
	        
	        int structured = result.getInt("structured_models_only");
	        if (structured == 1) {
	        	onlyStructured = true;
	        } else {
	        	onlyStructured = false;
	        }
	        result.close();
	        
	        loops = loadLoops(queryId);
		} catch (SQLException e) {
			System.out.println(sqlForQueryInfo);
			throw new QueryLoadingFailed("Could not load: " + queryId, e);
		}
		
        Query query = new ImmutableQuery(queryId, loops, onlyStructured, ilSet, hlSet, modelType);
        
        return query;
	}
	
	public void cleanUp() {
		try {
			connection.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}
}
