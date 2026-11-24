package morphmodels.evolution.likelihood;

import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

/**
 * Probability that a character is parsimony uninformative among a given set of scored taxa.
 *
 * A pattern is uninformative when one state is a background and every other state it uses appears
 * exactly once, so the sum runs over a background state and a set of singleton states. Enumerating
 * those patterns is exponential in the number of taxa, so they are summed by a recursion that
 * carries, for each state at a node, which singleton states the subtree below it has used.
 *
 * Unscored taxa take any state and contribute no singleton, which is what makes the result a
 * probability over the scored taxa alone.
 *
 * @author June Walker
 */
public class UninformativeProbability {

    /** mixture-weighted probability that the pattern is uninformative among the scored taxa */
    public static double uninformative(TreeInterface tree, SubstitutionModel sm,
                                       SiteModel.Base siteModel, BranchRateModel.Base clock,
                                       boolean[] scored, int stateCount) {
        return uninformative(tree, sm, siteModel, clock, scored, stateCount, false);
    }

    /**
     * With ordered = true the uninformative class is the one of the ordered (Wagner) model, where
     * a pattern is informative only when some threshold has two taxa on each side. Such a pattern
     * is uninformative exactly when one state is a background, at most one taxon lies strictly
     * below it and at most one strictly above, so the recursion carries two bits rather than a set.
     */
    public static double uninformative(TreeInterface tree, SubstitutionModel sm,
                                       SiteModel.Base siteModel, BranchRateModel.Base clock,
                                       boolean[] scored, int stateCount, boolean ordered) {
        int scoredCount = 0;
        for (boolean b : scored)
            if (b) scoredCount++;
        if (scoredCount < 1)
            return 1.0;
        // either criterion needs two taxa on each side, so fewer than four can never be informative
        if (scoredCount < 4)
            return 1.0;

        double total = 0.0;
        int categories = siteModel.getCategoryCount();
        double[] proportions = siteModel.getCategoryProportions(tree.getRoot());
        for (int category = 0; category < categories; category++) {
            double rate = siteModel.getRateForCategory(category, tree.getRoot());
            total += proportions[category]
                   * forRate(tree, sm, clock, scored, stateCount, scoredCount, rate, ordered);
        }
        return total;
    }

    private static double forRate(TreeInterface tree, SubstitutionModel sm, BranchRateModel.Base clock,
                                  boolean[] scored, int k, int scoredCount, double rate, boolean ordered) {
        int nodeCount = tree.getNodeCount();
        int masks = ordered ? 4 : 1 << (k - 1);
        double[] freqs = sm.getFrequencies();
        double[][] probs = new double[nodeCount][k * k];
        for (Node node : tree.getNodesAsArray()) {
            if (node.isRoot())
                continue;
            double branchRate = rate * (clock == null ? 1.0 : clock.getRateForBranch(node));
            sm.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(),
                    branchRate, probs[node.getNr()]);
        }

        double total = 0.0;
        for (int background = 0; background < k; background++) {
            double[][] partials = new double[nodeCount][k * masks];
            fill(tree.getRoot(), partials, probs, scored, k, masks, background, ordered);
            double[] root = partials[tree.getRoot().getNr()];
            for (int mask = 0; mask < masks; mask++) {
                double weight = ordered ? 1.0 : weight(mask, scoredCount);
                if (weight == 0.0)
                    continue;
                double p = 0.0;
                for (int state = 0; state < k; state++)
                    p += freqs[state] * root[state * masks + mask];
                total += weight * p;
            }
        }
        return total;
    }

    /**
     * An all-distinct pattern is reached once for every state that could be called the background,
     * so it is counted scoredCount times and divided back down. More singletons than scored taxa
     * is impossible.
     */
    private static double weight(int mask, int scoredCount) {
        int bits = Integer.bitCount(mask);
        if (bits >= scoredCount)
            return 0.0;
        if (bits < scoredCount - 1)
            return 1.0;
        return 1.0 / scoredCount;
    }

    /** post-order fill of the (state, singleton set) partials */
    private static void fill(Node node, double[][] partials, double[][] probs, boolean[] scored,
                             int k, int masks, int background, boolean ordered) {
        double[] here = partials[node.getNr()];

        if (node.isLeaf()) {
            if (!scored[node.getNr()]) {
                // any state, no singleton used
                for (int state = 0; state < k; state++)
                    here[state * masks] = 1.0;
                return;
            }
            here[background * masks] = 1.0;
            for (int state = 0; state < k; state++) {
                if (state == background)
                    continue;
                int bit = ordered ? (state < background ? 0 : 1) : bitOf(state, background);
                here[state * masks + (1 << bit)] = 1.0;
            }
            return;
        }

        for (Node child : node.getChildren())
            fill(child, partials, probs, scored, k, masks, background, ordered);

        // the branch message of each child, then a disjoint split of the singleton set between them
        int children = node.getChildCount();
        double[][] msg = new double[children][k * masks];
        for (int c = 0; c < children; c++) {
            Node child = node.getChild(c);
            double[] below = partials[child.getNr()];
            double[] p = probs[child.getNr()];
            for (int from = 0; from < k; from++)
                for (int mask = 0; mask < masks; mask++) {
                    double sum = 0.0;
                    for (int to = 0; to < k; to++)
                        sum += p[from * k + to] * below[to * masks + mask];
                    msg[c][from * masks + mask] = sum;
                }
        }

        for (int state = 0; state < k; state++)
            for (int mask = 0; mask < masks; mask++)
                here[state * masks + mask] = combine(msg, 0, children, state, mask, masks);
    }

    /** sum over ways of splitting mask disjointly among the children from index c onwards */
    private static double combine(double[][] msg, int c, int children, int state, int mask, int masks) {
        if (c == children - 1)
            return msg[c][state * masks + mask];
        double total = 0.0;
        int sub = mask;
        while (true) {
            total += msg[c][state * masks + sub]
                   * combine(msg, c + 1, children, state, mask & ~sub, masks);
            if (sub == 0)
                break;
            sub = (sub - 1) & mask;
        }
        return total;
    }

    /** bit position of a non-background state, the states being indexed in order with the background left out */
    private static int bitOf(int state, int background) {
        return state < background ? state : state - 1;
    }
}
