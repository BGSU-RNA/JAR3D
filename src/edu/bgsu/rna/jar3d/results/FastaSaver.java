package edu.bgsu.rna.jar3d.results;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class FastaSaver extends AbstractResultsSaver {
	
	private BufferedWriter writer;
	
	public FastaSaver(String filename) throws IOException {
		this.writer = new BufferedWriter(new FileWriter(filename));
	}

	public void save(LoopResult results) throws SaveFailed {

		String queryId = results.queryId();
		String loopId = new Integer(results.loopId()).toString();
		
		for (SequenceResult result: results.sequenceResults()) {
			String motifId = result.motifId();
			String score = new Double(result.score()).toString();
			String percentile = new Double(result.percentile()).toString();
			String editDistance = new Integer(result.editDistance()).toString();
			String rotated = new Boolean(result.isRotated()).toString();
			String line = join(queryId, loopId, motifId, score, percentile, editDistance, rotated);
			try {
				writer.write(line);
				writer.newLine();
			} catch (IOException e) {
				throw new SaveFailed(e);
			}
		}
	}
	
	private String join(String ... data) {
		StringBuilder builder = new StringBuilder();
		for(int i = 0; i < data.length; i++) {
			if (data[i].equals("")) {
				data[i] = "NA";
			}
			builder.append(data[i]);
			if (i + 1 != data.length) {
				builder.append(",");
			}
		}
		return builder.toString();
	}

	public void markFailure(String queryId) throws SaveFailed {
		throw new RuntimeException("Not implemented yet");
	}

	public void cleanUp() throws SaveFailed {
		try {
			this.writer.close();
		} catch (IOException e) {
			throw new SaveFailed(e);
		}
	}

}
