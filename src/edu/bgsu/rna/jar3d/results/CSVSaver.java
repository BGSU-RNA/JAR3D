package edu.bgsu.rna.jar3d.results;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class CSVSaver extends AbstractResultsSaver {
	
	private BufferedWriter loopWriter;
	
	private BufferedWriter sequenceWriter;
	
	public CSVSaver(String loopFile, String sequenceFile) throws IOException {
		this.loopWriter = new BufferedWriter(new FileWriter(loopFile));
		this.sequenceWriter = new BufferedWriter(new FileWriter(sequenceFile));
	}
	
	private String format(double number) {
		return String.format("%.5g", number);
	}
	
	private void saveLoopResults(LoopResult results) throws SaveFailed {
		String queryId = results.queryId();
		String loopId = new Integer(results.loopId()).toString();
		String motifId = results.modelId();
		String meanScore = format(results.meanScore());
		String medianScore = format(results.medianScore());
		String meanPercentile = format(results.meanPercentile());
		String medianPercentile = format(results.medianPercentile());
		String meanEditDistance = format(results.meanEditDistance());
		String medianEditDistance = format(results.medianEditDistance());
		
		String line = join(queryId, loopId, motifId, meanScore, medianScore, 
				meanPercentile, medianPercentile, meanEditDistance, medianEditDistance);
		
		try {
			loopWriter.write(line);
			loopWriter.newLine();
		} catch (IOException e) {
			throw new SaveFailed(e);
		}
	}
	
	private void saveSequenceResults(int loopId, SequenceResult result) throws SaveFailed {
		String queryId = result.queryId();
		String loopIdString = new Integer(loopId).toString();
		String sequenceId = result.sequenceId();
		String motifId = result.motifId();
		String score = format(result.score());
		String percentile = format(result.percentile());
		String editDistance = new Integer(result.editDistance()).toString();
		String rotated = new Boolean(result.isRotated()).toString();
		String line = join(queryId, loopIdString, sequenceId, motifId, score, percentile, editDistance, rotated);
		try {
			sequenceWriter.write(line);
			sequenceWriter.newLine();
		} catch (IOException e) {
			throw new SaveFailed(e);
		}
	}

	public void save(LoopResult results) throws SaveFailed {
		saveLoopResults(results);
		for(SequenceResult sequenceResult: results.sequenceResults()) {
			saveSequenceResults(results.loopId(), sequenceResult);
		}		
	}
	
	private String join(String ... data) {
		StringBuilder builder = new StringBuilder();
		for(int i = 0; i < data.length; i++) {
			if (data[i] == null || data[i].equals("")) {
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
			loopWriter.close();
			sequenceWriter.close();
		} catch (IOException e) {
			throw new SaveFailed(e);
		}
	}

	public void writeHeader() throws SaveFailed {
		String sequenceLine = join("queryId", "loopId", "sequenceId", "motifId", 
				"score", "percentile", "editDistance", "rotated");
		String loopLine = join("queryId", "loopId", "motifId", "meanScore", "medianScore", 
				"meanPercentile", "medianPercentile", "meanEditDistance", "medianEditDistance");
		
		try {
			sequenceWriter.write(sequenceLine);
			sequenceWriter.newLine();
			loopWriter.write(loopLine);
			loopWriter.newLine();
		} catch(IOException e) {
			throw new SaveFailed(e);
		}
		
	}

}
