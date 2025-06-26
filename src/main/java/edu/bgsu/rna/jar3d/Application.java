package edu.bgsu.rna.jar3d;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import edu.bgsu.rna.jar3d.io.loaders.QueryLoader;
import edu.bgsu.rna.jar3d.io.loaders.QueryLoadingFailed;
import edu.bgsu.rna.jar3d.io.writers.ResultSaver;
import edu.bgsu.rna.jar3d.io.writers.SaveFailed;
import edu.bgsu.rna.jar3d.loop.Loop;
import edu.bgsu.rna.jar3d.loop.LoopType;
import edu.bgsu.rna.jar3d.query.Query;
import edu.bgsu.rna.jar3d.results.LoopResult;
import edu.bgsu.rna.jar3d.results.SequenceResult;

/**
 * An Application is simple way to load sequences and models, run the models over the sequences
 * and save the results.
 * This can only load from the file system currently and uses only bp models at the moment.
 */
public class Application {

	/** Default range limit for how far to look left and right in the sequence, mostly for alignments */
	// TODO 2013-11-05 CLZ The user needs a way to override the default range limit
	// 100 will be adequate for motifs, but not for large alignments
	public static final int DEFAULT_RANGE_LIMIT = 100;

	/** Default model type */
	// 2013-11-05 CLZ This should no longer be used.  No default.
	public static final String DEFAULT_MODEL_TYPE = "bp";

	/** Default model version */
	// 2013-11-05 CLZ This should no longer be used.  No default.
	public static final String DEFAULT_VERSION = "0.6";

	/** The query loader to use. */
	private final QueryLoader loader;

	/** The object to save results with. */
	private final ResultSaver saver;

	/** The model type to use. */
	private final String modelType;

	/** The model version to use. */
	private final String version;

	/** The range limit to scan. */
	private final int rangeLimit;

	/**
	 * Build a new Application.
	 *
	 * @param loader The query loader to use.
	 * @param saver The result saver to use.
	 * @param modelType The model type to use.
	 * @param version The version to use.
	 * @param rangeLimit The range limit to scan.
	 */
	public Application(QueryLoader loader, ResultSaver saver, String modelType, String version, int rangeLimit) {
		this.loader = loader;
		this.saver = saver;
		this.modelType = modelType;
		this.version = version;
		this.rangeLimit = rangeLimit;
	}

	/**
	 * Create a new Application. The modelType, version, and rangeLimit are set to the default ones.
	 *
	 * @param loader The query loader to use.
	 * @param saver The results saver to use.
	 */
	public Application(QueryLoader loader, ResultSaver saver) {
		this(loader, saver, DEFAULT_MODEL_TYPE, DEFAULT_VERSION, DEFAULT_RANGE_LIMIT);
	}

	/**
	 * Load then run a query and return the results.
	 *
	 * @param queryId Query id to load.
	 * @param base Base path to the models.
	 * @return The results.
	 * @throws QueryLoadingFailed
	 */
	public List<List<LoopResult>> runQuery(String queryId, String base) throws QueryLoadingFailed {
		Query query = loader.load(queryId);
		return runQuery(query, base);
	}

	/**
	 * Load then run a query and return the results.
	 *
	 * @param queryId Query id to load.
	 * @param ILbase Base path to the IL models.
	 * @param HLbase Base path to the HL models.
	 * @return The results.
	 * @throws QueryLoadingFailed
	 * @throws SaveFailed
	 */

	public void runQuery(String queryId, String ILbase, String HLbase) throws QueryLoadingFailed, SaveFailed {
		Query query = loader.load(queryId);
		this.runQuery(query, ILbase, HLbase);
	}

	/**
	 * Run a query and return the results.
	 *
	 * @param base The base bath to the models.
	 * @param query The query to run.
	 * @return The results.
	 */
	public List<List<LoopResult>> runQuery(Query query, String base) {
		List<List<LoopResult>> allResults = new ArrayList<List<LoopResult>>();
		//Load Motif Group information
		File f = new File(base);
		String folder = f.getParent();

		System.out.println("Looking for a list of motifs to use in "+base);
		System.out.println("Looking for motifs in "+folder);

		System.setProperty("user.dir", folder);

		Vector<String> modelNames = Sequence.getModelNames(base, modelType, false);
		HashMap<String,MotifGroup> groupData = webJAR3D.loadMotifGroups(base, modelType);


		for(Loop loop: query) {
			List<LoopResult> results = motifParse(modelNames, groupData, loop, true);
			allResults.add(results);
		}
		return allResults;
	}

	/**
	 * Run a query and return the results.
	 *
	 * @param ILbase The base path to the IL models.
	 * @param HLbase The base path to the HL models
	 * @param query The query to run.
	 * @return The results.
	 * @throws SaveFailed
	 */

	public void runQuery(Query query, String ILbase, String HLbase) throws SaveFailed {

		List<LoopResult> results;

		for(Loop loop: query) {
			String type = loop.getLoopType().getShortName();
			if (type.equals("HL")){
				//Load HL Motif Group information
				File hlf = new File(HLbase);
				String hlfolder = hlf.getParent();

				System.out.println("Looking for a list of motifs to use in "+HLbase);
				System.out.println("Looking for motifs in "+hlfolder);

				System.setProperty("user.dir", hlfolder);

				Vector<String> HLModelNames = Sequence.getModelNames(HLbase, modelType, false);
				HashMap<String,MotifGroup> HLGroupData = webJAR3D.loadMotifGroups(HLbase, modelType);
				results = motifParse(HLModelNames, HLGroupData, loop, false);
			} else if (type.equals("IL")){
				//Load IL Motif Group information
				File ilf = new File(ILbase);
				String ilfolder = ilf.getParent();

				System.out.println("Looking for a list of motifs to use in "+ILbase);
				System.out.println("Looking for motifs in "+ilfolder);

				System.setProperty("user.dir", ilfolder);

				Vector<String> ILModelNames = Sequence.getModelNames(ILbase, modelType, false);
				HashMap<String,MotifGroup> ILGroupData = webJAR3D.loadMotifGroups(ILbase, modelType);
				results = motifParse(ILModelNames, ILGroupData, loop, false);
			}
		} else if (type.startsWith("J")){
			// Load J3 or J4 ... Motif Group information
			// Modify the HL path to get the JX path
			String JLbase = HLbase.replace("/HL/", "/" + type + "/");

			File jlf = new File(JLbase);
			String jlfolder = jlf.getParent();

			System.out.println("Looking for a list of motifs to use in "+JLbase);
			System.out.println("Looking for motifs in "+jlfolder);

			System.setProperty("user.dir", jlfolder);

			Vector<String> JLModelNames = Sequence.getModelNames(JLbase, modelType, false);
			HashMap<String,MotifGroup> JLGroupData = webJAR3D.loadMotifGroups(JLbase, modelType);
			results = motifParse(JLModelNames, JLGroupData, loop, false);
		}
		saver.save(results);
		}
		saver.cleanUp();
	}

	/**
	 * Run a loop query and return the results.
	 *
	 * @param queryId Query to load.
	 * @param loopNumber Which loop to run on.
	 * @param folder Folder the model files are in.
	 * @param model Name of the model to load.
	 * @return The results.
	 * @throws QueryLoadingFailed
	 * @throws SaveFailed
	 */

	private void runQuery(String queryID, String loopNumber, String folder,
			String model) throws QueryLoadingFailed, SaveFailed {
		Query query = loader.load(queryID);
		this.runQuery(query, loopNumber, folder, model);
	}

	/**
	 * Run a loop query and return the results.
	 *
	 * @param queryId Query to load.
	 * @param loopNumber Which loop to run on.
	 * @param folder Folder the model files are in.
	 * @param model Name of the model to load.
	 * @return The results.
	 * @throws SaveFailed
	 */

	private void runQuery(Query query, String loopNumber, String folder,
			String model) throws SaveFailed {
		// Load Motif Group
		MotifGroup Group = new MotifGroup(folder, modelType, model);
		Loop loop = query.getLoops().get(Integer.parseInt(loopNumber));
		List<SequenceResult> results;

		results = Alignment.doSingleDBQuery(loop, Group, model, rangeLimit);

		saver.saveSequenceResults(results,query,loopNumber);
	}

	/**
	 * Run a loop against a single loop and return the results. This will only score internal loops, all other loop
	 * types will give empty results.
	 *
	 * @param modellist The path to a file telling what models to use
	 * @param loop The loop.
	 * @return Results of running the loop.
	 */
	private List<LoopResult> motifParse(Vector<String> modelNames, HashMap<String,MotifGroup> groupData, Loop loop, boolean saveSeqRes) {
		List<LoopResult> result = new ArrayList<LoopResult>();

		if (modelNames.size() == 0) {
			System.out.println("Found " + modelNames.size() + " model files");
		}

		System.out.println("Application.motifParse: " +loop.getLoopType());

		result = Alignment.doLoopDBQuery(loop, modelNames, groupData, rangeLimit, saveSeqRes);

		return result;
	}

	/**
	 * Save results using the ResultSaver.
	 *
	 * @param results The results to save.
	 * @throws SaveFailed
	 */
	public void saveResults(List<List<LoopResult>> results) throws SaveFailed {
		for(List<LoopResult> res: results) {
			saver.save(res);
		}
		saver.cleanUp();
	}

	/**
	 * Load a query, run the query and then save the results.
	 *
	 * @param queryId Query to load.
	 * @param base Base path to the models.
	 * @throws SaveFailed
	 * @throws QueryLoadingFailed
	 */

	public void runAndSave(String queryId, String base) throws SaveFailed, QueryLoadingFailed {
		saver.writeHeader();
		List<List<LoopResult>> results = this.runQuery(queryId, base);
		for(List<LoopResult> res: results) {
			saver.save(res);
		}
		saver.cleanUp();
	}

	/**
	 * Load a query, run the query and then save the results.
	 *
	 * @param queryId Query to load.
	 * @param ILbase Base path to the IL models.
	 * @param HLbase Base path to the HL models.
	 * @throws SaveFailed
	 * @throws QueryLoadingFailed
	 */

	public void runAndSave(String queryId, String ILbase, String HLbase) throws SaveFailed, QueryLoadingFailed {
		saver.writeHeader();
		this.runQuery(queryId, ILbase, HLbase);
	}

	/**
	 * Load a loop query, run the query and then save the results.
	 *
	 * @param queryId Query to load.
	 * @param loopNumber Which loop to run on.
	 * @param folder Folder the model files are in.
	 * @param model Name of the model to load.
	 * @throws SaveFailed
	 * @throws QueryLoadingFailed
	 */

	public void runAndSave(String queryId, String loopNumber, String folder, String model) throws SaveFailed, QueryLoadingFailed {
		saver.writeHeader();
		this.runQuery(queryId, loopNumber, folder, model);
	}

}
