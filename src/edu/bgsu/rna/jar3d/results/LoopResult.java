package edu.bgsu.rna.jar3d.results;

import java.util.List;

public interface LoopResult {

//	"insert into bygroup (id, meanscore, meanpercentile, meaneditdist, medianscore, medianpercentile, medianeditdist, signature, rotation, groupnum) values('%s', %f, %f, %f, %f, %f, %f,'%s',%d,%s) ", id,modelScores[g],meanQuant,meanMinDist,medianLL,medianQuant,medianMinDist,sig,reversed[index],groupName);

	public int loopId();

	public String modelId();

	public double meanEditDistance();

	public double meanScore();

	public double meanPercentile();

	public double medianScore();

	public double medianPercentile();

	public double medianEditDistance();

	public String signature();

	public boolean isRotated();

	public List<SequenceResult> sequenceResults();
}
