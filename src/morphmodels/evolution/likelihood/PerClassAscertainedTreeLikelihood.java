package morphmodels.evolution.likelihood;

import beast.base.core.Description;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.tree.Node;

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

    protected PerClassAscertainedAlignment perClass;
    protected Alignment data;
    private double[] byClass, byPattern;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        data = dataInput.get();
        if (useAscertainedSitePatterns && data instanceof PerClassAscertainedAlignment)
            perClass = (PerClassAscertainedAlignment) data;
    }

    @Override
    protected void calcLogP() {
        if (perClass == null) {
            super.calcLogP();
            return;
        }
        logP = 0.0;
        final double[] correction = perClass.isInformative()
                ? informativeCorrections() : perClass.getAscertainmentCorrections(patternLogLikelihoods);
        for (int i = 0; i < data.getPatternCount(); i++)
            logP += (patternLogLikelihoods[i] - correction[i]) * data.getPatternWeight(i);
    }

    /**
     * log(1 - P(uninformative)) for each pattern's class. There are too many uninformative
     * patterns to append and read back, so the probability is summed by recursion over the tree.
     */
    protected double[] informativeCorrections() {
        java.util.List<boolean[]> classes = perClass.getScoredOfClass();
        int[] classOf = perClass.getClassOfPattern();
        int states = data.getDataType().getStateCount();
        if (byClass == null || byClass.length != classes.size())
            byClass = new double[classes.size()];
        if (byPattern == null || byPattern.length != classOf.length)
            byPattern = new double[classOf.length];

        for (int cl = 0; cl < classes.size(); cl++) {
            boolean[] scoredByNode = new boolean[treeInput.get().getNodeCount()];
            for (Node node : treeInput.get().getNodesAsArray())
                if (node.isLeaf()) {
                    int taxon = data.getTaxonIndex(node.getID());
                    scoredByNode[node.getNr()] = taxon >= 0 && classes.get(cl)[taxon];
                }
            double unobserved = UninformativeProbability.uninformative(treeInput.get(),
                    substitutionModel, m_siteModel, branchRateModel, scoredByNode, states);
            byClass[cl] = Math.log(1.0 - Math.min(unobserved, 1.0 - 1e-15));
        }
        for (int pattern = 0; pattern < classOf.length; pattern++)
            byPattern[pattern] = classOf[pattern] < 0 ? 0.0 : byClass[classOf[pattern]];
        return byPattern;
    }
}
