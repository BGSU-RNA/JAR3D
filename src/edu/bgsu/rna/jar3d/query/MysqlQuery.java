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
public class MysqlQuery extends AbstractQuery {
    private Connection connection;

    private PreparedStatement sqlForQueryInfo;
    
    private PreparedStatement sqlForLoops;
        
    private List<Loop> loops = new ArrayList<Loop>();
    
    private boolean onlyStructured = false;
    
    private String ilSet = null;
    
    private String hlSet = null;
    
    private String modelType = null;
    
    private String id;
    
    public MysqlQuery(String username, String password, String dbConnection, String queryId) throws SQLException {
        connection = DriverManager.getConnection(dbConnection, username, password);
        String querySql = "SELECT group_set, model_type, structured_models_only FROM query_info WHERE query_id = ?;";
        String loopSql = "SELECT id, loop_sequence, loop_type FROM query_sequences WHERE query_id = ? and loop_id = ?;";
        sqlForQueryInfo = connection.prepareStatement(querySql);
        sqlForLoops = connection.prepareStatement(loopSql);
        id = queryId;
        loadData();
        connection.close();
    }
    
    private Loop loadLoop(int index) throws SQLException {
    	sqlForLoops.setString(1, id);
    	sqlForLoops.setLong(2, index);
    	ResultSet results = sqlForLoops.executeQuery();
    	long id = -1;
    	
    	String type = null;
    	List<String> sequences = new ArrayList<String>();
    	while (results.next()) {
    		String sequence = results.getString("loop_sequence");
    		type = results.getString("loop_type");
    		id = results.getLong("id");
    		sequences.add(sequence);
    	}
    	results.close();

    	return new BasicLoop(id, sequences, type);
    }

    private void loadLoops() throws SQLException {
        String loopCountSql = "SELECT MAX(loop_id) AS max FROM query_sequences where query_id = ?;";
        PreparedStatement sqlForLoopCount = connection.prepareStatement(loopCountSql);
        sqlForLoopCount.setString(1, id);
        
        ResultSet result = sqlForLoopCount.executeQuery();
        boolean found = result.first();
        if (!found) {
        	throw new IndexOutOfBoundsException("Could not find any loops for query: " + id);
        }
        
        int loopCount = result.getInt("max") + 1;
        result.close();
        
        for (int i = 0; i < loopCount; i++) {
        	loops.add(loadLoop(i));
        }   
    }
    
    private void loadMetaData() throws SQLException {
        sqlForQueryInfo.setString(1, id);
        
        ResultSet result = sqlForQueryInfo.executeQuery();
        boolean found = result.first();
        
        if (!found) {
        	throw new IndexOutOfBoundsException("Could not find query with id: " + id);
        }
        
        modelType = result.getString("model_type");
        String groupSet = result.getString("group_set");
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
    }
    
    private void loadData() throws SQLException {
        loadMetaData();
        loadLoops();
    }

	public String modelType() {
		return modelType;
	}

	public boolean onlyStructured() {
		return onlyStructured;
	}

	public List<Loop> getLoops() {
		return loops;
	}

	public String getId() {
		return id;
	}

	public String getILSetName() {
		return ilSet;
	}

	public String getHLSetName() {
		return hlSet;
	}
}
