package edu.bgsu.rna.jar3d.io.writers;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import edu.bgsu.rna.jar3d.results.LoopResult;
import edu.bgsu.rna.jar3d.results.SequenceResult;
import edu.bgsu.rna.jar3d.query.Query;

/**
 * This class saves results to a CSV file.
 */
public class CSVSaver extends AbstractResultsSaver {

    /** The Writer for the aggregated loop data. */
	private BufferedWriter loopWriter;

    /** The Writer for the sequence level data. */
	private BufferedWriter sequenceWriter;

	public CSVSaver(String loopFile, String sequenceFile) throws IOException {
		this.loopWriter = new BufferedWriter(new FileWriter(loopFile));
		this.sequenceWriter = new BufferedWriter(new FileWriter(sequenceFile));
	}

	private String format(double number) {
		return String.format("%.5g", Math.max(-9999,number));
	}

	private void saveLoopResults(LoopResult results) throws SaveFailed {
		//need to strip commas so file can be read as a csv
		String loopId = results.getLoop().getName().replace(",",";");
		String motifId = results.modelId();
		String meanScore = format(results.meanScore());
		String medianScore = format(results.medianScore());
		String meanInteriorEditDistance = format(results.meanInteriorEditDistance());
		String medianInteriorEditDistance = format(results.medianInteriorEditDistance());
		String meanFullEditDistance = format(results.meanFullEditDistance());
		String medianFullEditDistance = format(results.medianFullEditDistance());
		String rotation = format(results.bestRotation());
		String meanCutoff = format(results.meanCutoff());
		String meanCutoffScore = format(results.meanCutoffScore());
		
		String line = join(loopId, motifId, meanCutoff, meanCutoffScore, meanScore, medianScore,
				meanInteriorEditDistance, medianInteriorEditDistance, meanFullEditDistance, medianFullEditDistance, rotation);

		try {
			loopWriter.write(line);
			loopWriter.newLine();
		} catch (IOException e) {
			throw new SaveFailed(e);
		}
	}

	private void saveSequenceResults(String loopName, SequenceResult result) throws SaveFailed {
		//need to replace commas so file can be read as a csv
		String sequenceId = String.valueOf(result.sequenceId());
		String motifId = result.motifId();
		String score = format(result.score());
		String interiorEditDistance = Integer.valueOf(result.InteriorEditDistance()).toString();
		String fullEditDistance = Integer.valueOf(result.FullEditDistance()).toString();
		String rotated = format(result.bestRotation());
		String cutoff = String.valueOf(result.cutoff());
		String cutoffscore = String.valueOf(result.cutoffscore());
		String line = join(loopName, sequenceId, motifId, cutoff, cutoffscore, score, interiorEditDistance, fullEditDistance, rotated);
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
			saveSequenceResults(results.getLoop().getName(), sequenceResult);
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
	
	// Does nothing on purpose
	public void markAllDone(String queryId) throws SaveFailed { }


	public void writeHeader() throws SaveFailed {
		String sequenceLine = join("filename", "sequenceId", "motifId", "passedCutoff", "meanCutoffScore", 
				"score", "interiorEditDistance", "fullEditDistance", "rotation");
		String loopLine = join("filename", "motifId", "%passedCutoff", "meanCutoffScore", "meanScore", "medianScore",
				"meanInteriorEditDistance", "medianInteriorEditDistance", "meanFullEditDistance", 
				"medianFullEditDistance", "rotation");
		
		try {
			sequenceWriter.write(sequenceLine);
			sequenceWriter.newLine();
			loopWriter.write(loopLine);
			loopWriter.newLine();
		} catch(IOException e) {
			throw new SaveFailed(e);
		}

	}

	public void saveSequenceResults(List<SequenceResult> results, Query query)
			throws SaveFailed {
		// TODO Auto-generated method stub
		
	}

}
