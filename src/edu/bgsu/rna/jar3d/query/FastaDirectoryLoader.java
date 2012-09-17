package edu.bgsu.rna.jar3d.query;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FastaDirectoryLoader implements QueryLoader {

    private FilenameFilter fastaFilter = new FilenameFilter() {
        public boolean accept(File dir, String name) {
            String lowerName = name.toLowerCase();
            return lowerName.endsWith(".fasta");
        }
    };
    
    private List<File> files;
    
    public FastaDirectoryLoader(File file) throws QueryLoadingFailed {
        if (!file.isDirectory()) {
            throw new QueryLoadingFailed("Must give a directory");
        }
        File[] files = file.listFiles(fastaFilter);
        this.files = Arrays.asList(files);
    }

    public Query load(String queryId) throws QueryLoadingFailed {
        List<Loop> loops = new ArrayList<Loop>();
        for (File file: files) {
            try {
                FastaLoader loader = new FastaLoader(file.getAbsolutePath());
                Query query = loader.load(null);
                loops.addAll(query.getLoops());
                loader.cleanUp();
            } catch (IOException e) {
                throw new QueryLoadingFailed(e);
            }
        }
        return new ImmutableQuery(queryId, loops, true, "IL", "HL", "bp_models");
    }

    public void cleanUp() throws QueryLoadingFailed { }
}
