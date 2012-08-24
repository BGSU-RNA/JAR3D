package edu.bgsu.rna.jar3d.query;

public interface QueryLoader {

	/**
	 * Load a query. 
	 * 
	 * @param queryId The id of the query to load.
	 * @return The requested query. 
	 */
	public Query load(String queryId) throws QueryLoadingFailed;

	/**
	 * Perform any clean up needed, such as closing connections etc.
	 */
	public void cleanUp() throws QueryLoadingFailed;
}
