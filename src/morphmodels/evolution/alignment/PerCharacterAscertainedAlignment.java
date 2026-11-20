package morphmodels.evolution.alignment;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Groups patterns by their sets of scored taxa after the usual pattern calculation.
 * Appends one constant pattern per state per class, or one per class with exchangeableStates,
 * leaving unscored taxa as missing. Appended patterns have zero weight and do not count as data.
 * Set only ascertained="true" and use PerCharacterAscertainedTreeLikelihood.
 * If most characters have distinct sets of scored taxa, this multiplies the pattern count.
 *
 * @author June Walker
 */
@Description("Provides ascertainment corrections conditioned on variability within each pattern's "
        + "set of scored taxa")
public class PerCharacterAscertainedAlignment extends Alignment {

    final public Input<Boolean> exchangeableStatesInput = new Input<>("exchangeableStates",
            "if true, append one constant pattern per class instead of one per state, and multiply its "
            + "probability by the number of states. Valid when the states are exchangeable, as under "
            + "LewisMK with equal frequencies, and not under an ordered model", false);

    /** class index for each pattern, or -1 for an appended constant pattern */
    protected int[] classOfPattern;
    /** indices of the appended constant patterns of each class */
    protected List<int[]> constantPatternsOfClass;
    /** reused between calls so no array is allocated per likelihood evaluation */
    protected double[] correctionOfClass, correctionOfPattern;
    /** number of constant patterns each appended one stands for, 1 unless states are exchangeable */
    protected int stateMultiplier = 1;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        if (isAscertainedInput.get())
            buildClasses();
    }

    /**
     * group patterns by their sets of scored taxa and append zero-weight constant patterns for each class
     */
    protected void buildClasses() {
        final int taxonCount = taxaNames.size();
        final int stateCount = m_dataType.getStateCount();
        if (stateCount < 2)
            throw new IllegalArgumentException("expected a data type with a known number of states, got "
                    + stateCount);
        final int missingCode = missingCode(stateCount);
        final int patternCount = getPatternCount();

        Map<String, Integer> classOfKey = new HashMap<>();
        List<boolean[]> scoredOfClass = new ArrayList<>();
        classOfPattern = new int[patternCount];

        for (int pattern = 0; pattern < patternCount; pattern++) {
            boolean[] scored = new boolean[taxonCount];
            StringBuilder key = new StringBuilder(taxonCount);
            for (int taxon = 0; taxon < taxonCount; taxon++) {
                scored[taxon] = !isUnscored(sitePatterns[pattern][taxon], stateCount);
                key.append(scored[taxon] ? '1' : '0');
            }
            Integer cl = classOfKey.get(key.toString());
            if (cl == null) {
                cl = scoredOfClass.size();
                classOfKey.put(key.toString(), cl);
                scoredOfClass.add(scored);
            }
            classOfPattern[pattern] = cl;
        }

        final boolean exchangeable = exchangeableStatesInput.get();
        final int perClass = exchangeable ? 1 : stateCount;
        final int appended = scoredOfClass.size() * perClass;
        int[][] patterns = Arrays.copyOf(sitePatterns, patternCount + appended);
        int[] weights = Arrays.copyOf(patternWeight, patternCount + appended);
        int[] classes = Arrays.copyOf(classOfPattern, patternCount + appended);
        constantPatternsOfClass = new ArrayList<>();

        int next = patternCount;
        for (boolean[] scored : scoredOfClass) {
            int[] mine = new int[perClass];
            for (int state = 0; state < perClass; state++) {
                int[] pattern = new int[taxonCount];
                for (int taxon = 0; taxon < taxonCount; taxon++)
                    pattern[taxon] = scored[taxon] ? state : missingCode;
                patterns[next] = pattern;
                weights[next] = 0;
                classes[next] = -1;
                mine[state] = next++;
            }
            constantPatternsOfClass.add(mine);
        }
        sitePatterns = patterns;
        patternWeight = weights;
        classOfPattern = classes;
        correctionOfClass = new double[scoredOfClass.size()];
        correctionOfPattern = new double[classes.length];
        stateMultiplier = exchangeable ? stateCount : 1;

        Log.info.println((getID() == null ? "alignment" : getID()) + ": " + patternCount + " patterns in "
                + scoredOfClass.size() + " classes of scored taxa, with " + appended
                + " constant patterns appended");
    }

    /** true if the code leaves the state unknown, so the taxon was not scored */
    protected boolean isUnscored(int code, int stateCount) {
        int[] states = m_dataType.getStatesForCode(code);
        return states == null || states.length >= stateCount;
    }

    /** find the missing code by asking the data type which code maps to every state */
    protected int missingCode(int stateCount) {
        for (int code = 0; ; code++) {
            int[] states;
            try {
                states = m_dataType.getStatesForCode(code);
            } catch (RuntimeException e) {
                break;
            }
            if (states == null)
                break;
            if (states.length >= stateCount)
                return code;
        }
        throw new IllegalArgumentException("no code for an unscored cell in data type "
                + m_dataType.getTypeDescription());
    }

    /**
     * return log(1 - summed probability of the class's constant patterns) for each pattern,
     * and zero for the appended constant patterns
     */
    public double[] getAscertainmentCorrections(double[] patternLogProbs) {
        for (int cl = 0; cl < correctionOfClass.length; cl++) {
            double unobserved = 0.0;
            for (int pattern : constantPatternsOfClass.get(cl))
                unobserved += Math.exp(patternLogProbs[pattern]);
            unobserved *= stateMultiplier;
            correctionOfClass[cl] = Math.log(1.0 - Math.min(unobserved, 1.0 - 1e-15));
        }
        for (int pattern = 0; pattern < correctionOfPattern.length; pattern++)
            correctionOfPattern[pattern] = classOfPattern[pattern] < 0
                    ? 0.0 : correctionOfClass[classOfPattern[pattern]];
        return correctionOfPattern;
    }

    public int getClassCount() {
        return constantPatternsOfClass == null ? 0 : constantPatternsOfClass.size();
    }
}
