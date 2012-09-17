package edu.bgsu.rna.jar3d.query;

import java.io.File;
import java.io.IOException;

public class GenericFastaLoader implements QueryLoader {

    private QueryLoader loader;
    
    public GenericFastaLoader(String filename) throws IOException, QueryLoadingFailed {
        File file = new File(filename);
        if (file.isFile()) {
            loader = new FastaLoader(filename);
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
