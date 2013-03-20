package edu.bgsu.rna.jar3d.io.loaders;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.bgsu.rna.jar3d.Sequence;
import edu.bgsu.rna.jar3d.loop.BasicLoop;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.ImmutableQuery;
import edu.bgsu.rna.jar3d.query.Query;

/**
 * This class loads a single loop from a FASTA formatted file. This assumes that the file represents a single loop. 
 */
public class FastaFileLoader implements QueryLoader {

	/** Reader for the file. */
	private BufferedReader reader;
	
	/** Filename we are reading from. */
	private final String filename;

	/**
	 * Create a new FastaFileLoader.
	 * 
	 * @param filename File to read from.
	 * @throws IOException
	 */
	public FastaFileLoader(String filename) throws IOException {
		this.filename = filename;
		this.reader = new BufferedReader(new FileReader(filename));
	}

	/**
	 * Load a Query from this FastaFileLoader. The given query id is used only as the query id of the returned query.
	 */
	public Query load(String queryId) throws QueryLoadingFailed {
		List<Sequence> sequences = new ArrayList<Sequence>();
		try {
			String line;
			String header = null;
			String sequence = null;

			while ((line = reader.readLine()) != null) {
				if (line.charAt(0) == '>') {
					if (header != null) {
						sequences.add(new Sequence(header, sequence));
					}
					header = line.substring(1);
				} else {
					sequence = line;
				}
			}
			sequences.add(new Sequence(header, sequence));
		} catch(IOException e) {
			throw new QueryLoadingFailed(e);
		}

		List<Loop> loops = new ArrayList<Loop>();
		loops.add(new BasicLoop(filename, 1, sequences));

		Query query = new ImmutableQuery(queryId, loops, true, "IL0.8", "", "IL");
		for(Loop loop: query) {
			loop.setQuery(query);
		}
		return query;
	}

	public void cleanUp() throws QueryLoadingFailed {
		try {
			reader.close();
		} catch (IOException e) {
			throw new QueryLoadingFailed(e);
		}
	}

}
