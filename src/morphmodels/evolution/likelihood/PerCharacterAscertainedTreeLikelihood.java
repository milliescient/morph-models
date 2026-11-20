package morphmodels.evolution.likelihood;

import beast.base.core.Description;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;

import morphmodels.evolution.alignment.PerCharacterAscertainedAlignment;

/**
 * Overrides calcLogP to subtract the ascertainment correction for each pattern's class.
 * Behaves as TreeLikelihood when the alignment is not a PerCharacterAscertainedAlignment.
 *
 * @author June Walker
 */
@Description("Applies ascertainment corrections by scored-taxon class for a "
        + "PerCharacterAscertainedAlignment, otherwise behaving as TreeLikelihood")
public class PerCharacterAscertainedTreeLikelihood extends TreeLikelihood {

    protected PerCharacterAscertainedAlignment perCharacter;
    protected Alignment data;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        data = dataInput.get();
        if (useAscertainedSitePatterns && data instanceof PerCharacterAscertainedAlignment)
            perCharacter = (PerCharacterAscertainedAlignment) data;
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
