
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;

import java.util.Map;
import java.util.HashMap;

/**
 * This is a class to load query data from the database.
 */
public class MysqlLoader {
    private Connection connection;

    private PreparedStatement sqlQuery;

    public MysqlLoader(String username, String password, String dbConnection) {
        connection = DriverManager.getConnection(dbConnection, username, password);
        String sql = "SELECT * from queries where query_id = ? limit 1";
        sqlQuery = connection.prepareStatement(sql);
    }

    public Map<String, String> getQuery(String queryId) {
        Map<String, String> query = new HashMap<String, String>();

        sqlQuery.setString(1, queryId);

        ResultSet result = sqlQuery.executeQuery();
        ResultSetMetaData meta = result.getMetaData();
        result.first();

        for (int i = 0; i < meta.getColumnCount(); i++) {
            String name = meta.getColumnName(i);
            String value = result.getString(i);
            query.put(name, value);
        }

        result.close();
        result.clearParameters();

        return query;
    }

    public close() {
      connection.close();
    }
}
