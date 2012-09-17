package edu.bgsu.rna.jar3d.query;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class FastaLoader implements QueryLoader {

	
	private BufferedReader reader;
	
	public FastaLoader(String filename) throws IOException {
		this.reader = new BufferedReader(new FileReader(filename));
	}
	
	public Query load(String queryId) throws QueryLoadingFailed {
    	List<String> sequences = new ArrayList<String>();
    	String name;
		try {
			String line;
			String header = null;
			String sequence = null;
			
			while ((line = reader.readLine()) != null) {
				if (line.charAt(0) == '>') {
					if (header != null) {
						sequences.add(sequence);
					}
					header = line.substring(1);
				} else {
					sequence = line;
				}
			}
			sequences.add(sequence);
			name = header;
		} catch(IOException e) {
			throw new QueryLoadingFailed(e);
		}
		
		List<Loop> loops = new ArrayList<Loop>();
		loops.add(new BasicLoop(name, 1, sequences, "IL"));
		
        return new ImmutableQuery(queryId, loops, true, "IL0.8", "", "IL");
	}

	public void cleanUp() throws QueryLoadingFailed {
		try {
			reader.close();
		} catch (IOException e) {
			throw new QueryLoadingFailed(e);
		}
	}

}
