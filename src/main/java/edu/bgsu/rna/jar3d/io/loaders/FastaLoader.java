package edu.bgsu.rna.jar3d.io.loaders;

import java.io.File;
import java.io.IOException;

import edu.bgsu.rna.jar3d.query.Query;

public class FastaLoader implements QueryLoader {

    private QueryLoader loader;
    
    public FastaLoader(String filename) throws IOException, QueryLoadingFailed {
        File file = new File(filename);
        if (file.isFile()) {
            loader = new FastaFileLoader(filename);
        } else {
            loader = new FastaDirectoryLoader(file);
        }
    }
    
    public Query load(String queryId) throws QueryLoadingFailed {
        return loader.load(queryId);
    }

    public void cleanUp() throws QueryLoadingFailed {
        loader.cleanUp();
    }

}
