package edu.bgsu.rna.jar3d.io.loaders;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.query.ImmutableQuery;
import edu.bgsu.rna.jar3d.query.Query;

/**
 * This class creates a query with from all fasta files in a directory. Each file is a separate loop but all are part
 * of the same query. 
 */
public class FastaDirectoryLoader implements QueryLoader {

	/** Should we load only structure models. */
    private boolean structuredOnly = true;

    /** The files we are going to load from. */
    private List<File> files;

    /** Filter to ensure we only read files that end with .fasta */
    private FilenameFilter fastaFilter = new FilenameFilter() {
        public boolean accept(File dir, String name) {
            String lowerName = name.toLowerCase();
            return lowerName.endsWith(".fasta");
        }
    };

    /**
     * Create a FastaDirectoryLoader.
     * 
     * @param file Path to the directory to read from.
     * @throws QueryLoadingFailed
     */
    public FastaDirectoryLoader(File file) throws QueryLoadingFailed {
        if (!file.isDirectory()) {
            throw new QueryLoadingFailed("Must give a directory");
        }
        File[] files = file.listFiles(fastaFilter);
        this.files = Arrays.asList(files);
    }

    /**
     * Load the query. The query id is used only to set the query id of the returned query.
     */
    public Query load(String queryId) throws QueryLoadingFailed {

        List<Loop> loops = new ArrayList<Loop>();
        for (File file: files) {
            try {
                FastaFileLoader loader = new FastaFileLoader(file.getAbsolutePath());
                Query query = loader.load(queryId);
                loops.addAll(query.getLoops());
                loader.cleanUp();
            } catch (IOException e) {
                throw new QueryLoadingFailed(e);
            }
        }

        Query query = new ImmutableQuery(queryId, loops, structuredOnly, "IL0.8", "", "IL");
        for(Loop loop: query) {
        	loop.setQuery(query);
        }
        return query;
    }

    /**
     * Does nothing.
     */
    public void cleanUp() throws QueryLoadingFailed { }
}
