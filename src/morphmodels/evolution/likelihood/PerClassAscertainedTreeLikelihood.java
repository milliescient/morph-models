package morphmodels.evolution.likelihood;

import beast.base.core.Description;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;

import morphmodels.evolution.alignment.PerClassAscertainedAlignment;

/**
 * Overrides calcLogP to subtract the ascertainment correction for each pattern's class.
 * Behaves as TreeLikelihood when the alignment is not a PerClassAscertainedAlignment.
 *
 * @author June Walker
 */
@Description("Applies ascertainment corrections by scored-taxon class for a "
        + "PerClassAscertainedAlignment, otherwise behaving as TreeLikelihood")
public class PerClassAscertainedTreeLikelihood extends TreeLikelihood {

    protected PerClassAscertainedAlignment perCharacter;
    protected Alignment data;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        data = dataInput.get();
        if (useAscertainedSitePatterns && data instanceof PerClassAscertainedAlignment)
            perCharacter = (PerClassAscertainedAlignment) data;
    }

    @Override
    protected void calcLogP() {
        if (perCharacter == null) {
            super.calcLogP();
            return;
        }
        logP = 0.0;
        final double[] correction = perCharacter.getAscertainmentCorrections(patternLogLikelihoods);
        for (int i = 0; i < data.getPatternCount(); i++)
            logP += (patternLogLikelihoods[i] - correction[i]) * data.getPatternWeight(i);
    }
}
