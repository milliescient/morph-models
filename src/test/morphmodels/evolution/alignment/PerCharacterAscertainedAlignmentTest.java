package test.morphmodels.evolution.alignment;

import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.StandardData;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;

import morphmodels.evolution.alignment.PerCharacterAscertainedAlignment;
import morphmodels.evolution.likelihood.PerCharacterAscertainedTreeLikelihood;
import morphmodels.evolution.substitutionmodel.LewisMK;
import morphmodels.evolution.substitutionmodel.Ordinal;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

/**
 * Six taxa on a fixed tree, ten characters with eight distinct sets of scored taxa.
 * The expected values were computed independently in Python and in MrBayes, which agree
 * with each other to the precision asserted here.
 *
 * @author June Walker
 */
public class PerCharacterAscertainedAlignmentTest {

    static final double PRECISION = 1e-5;
    static final String[] TAXA = {"A", "B", "C", "D", "E", "F"};
    static final String NEWICK = "(A:0.1,B:0.2,((C:0.15,D:0.05):0.12,(E:0.3,F:0.08):0.2):0.07);";

    static final String[] BINARY = {"00000?0011", "11?11?1?1?", "0?00000000",
                                    "111?111100", "000000??11", "1111?1?100"};
    static final String[] THREE = {"00202?0121", "2?111?2?22", "10?202222?",
                                   "122?001020", "0222?10?00", "2101?1?012"};
    static final String[] FOUR = {"30030?1302", "2?222?1?13", "21?131203?",
                                  "323?133130", "0330?00?20", "1013?2?211"};
    static final String[] THREE_ORDERED = {"01201?1002", "2?122?2?11", "21?100011?",
                                           "010?110020", "0000?20?20", "1200?0?210"};
    static final String[] FOUR_ORDERED = {"21330?1023", "2?002?3?12", "32?112223?",
                                          "101?300320", "1302?32?01", "0221?1?111"};

    @Test
    public void testBinary() throws Exception {
        assertEquals(-45.857762, perCharacter(BINARY, 2, false, false), PRECISION);
    }

    @Test
    public void testThreeStates() throws Exception {
        assertEquals(-69.586726, perCharacter(THREE, 3, false, false), PRECISION);
    }

    @Test
    public void testFourStates() throws Exception {
        assertEquals(-95.453812, perCharacter(FOUR, 4, false, false), PRECISION);
    }

    @Test
    public void testThreeStatesOrdered() throws Exception {
        assertEquals(-73.468901, perCharacter(THREE_ORDERED, 3, true, false), PRECISION);
    }

    @Test
    public void testFourStatesOrdered() throws Exception {
        assertEquals(-115.125951, perCharacter(FOUR_ORDERED, 4, true, false), PRECISION);
    }

    /** the correction of a class with exchangeable states needs only one of its constant patterns */
    @Test
    public void testExchangeableStatesAgree() throws Exception {
        for (int k = 2; k <= 4; k++) {
            String[] rows = k == 2 ? BINARY : k == 3 ? THREE : FOUR;
            assertEquals(perCharacter(rows, k, false, false), perCharacter(rows, k, false, true), PRECISION);
        }
    }

    /** an ordered model has an end state and an interior state, so its constant patterns differ */
    @Test
    public void testExchangeableStatesRejectOrdered() throws Exception {
        assertNotEquals(perCharacter(FOUR_ORDERED, 4, true, false),
                        perCharacter(FOUR_ORDERED, 4, true, true), 1e-3);
    }

    /** one class per distinct set of scored taxa, and one appended pattern per state in each */
    @Test
    public void testClasses() throws Exception {
        PerCharacterAscertainedAlignment aln = alignment(BINARY, 2, false);
        assertEquals(8, aln.getClassCount());
        assertEquals(10 + 8 * 2, aln.getPatternCount());
    }

    /** conditioning on fully scored constant patterns, as an Alignment does, gives a different answer */
    @Test
    public void testFullyScoredDiffers() throws Exception {
        List<Sequence> seqs = new ArrayList<>();
        for (int i = 0; i < TAXA.length; i++) {
            Sequence s = new Sequence();
            s.initByName("taxon", TAXA[i], "value", BINARY[i] + "01");
            seqs.add(s);
        }
        StandardData dt = new StandardData();
        dt.initByName("nrOfStates", 2);
        Alignment aln = new Alignment();
        aln.initByName("sequence", seqs, "dataType", "standard", "userDataType", dt,
                "ascertained", true, "excludefrom", 10, "excludeto", 12);
        assertEquals(-47.159019, logP(new TreeLikelihood(), aln, dt, false), PRECISION);
    }

    /** on a plain Alignment there are no classes to use, so the likelihood is a TreeLikelihood */
    @Test
    public void testPlainAlignmentUnchanged() throws Exception {
        StandardData dt = new StandardData();
        dt.initByName("nrOfStates", 2);
        List<Sequence> seqs = new ArrayList<>();
        for (int i = 0; i < TAXA.length; i++) {
            Sequence s = new Sequence();
            s.initByName("taxon", TAXA[i], "value", BINARY[i]);
            seqs.add(s);
        }
        Alignment plain = new Alignment();
        plain.initByName("sequence", seqs, "dataType", "standard", "userDataType", dt);
        assertEquals(logP(new TreeLikelihood(), plain, dt, false),
                     logP(new PerCharacterAscertainedTreeLikelihood(), plain, dt, false), PRECISION);
    }

    double perCharacter(String[] rows, int nStates, boolean ordered, boolean exchangeable) throws Exception {
        StandardData dt = new StandardData();
        dt.initByName("nrOfStates", nStates);
        return logP(new PerCharacterAscertainedTreeLikelihood(), alignment(rows, nStates, exchangeable), dt, ordered);
    }

    PerCharacterAscertainedAlignment alignment(String[] rows, int nStates, boolean exchangeable) throws Exception {
        List<Sequence> seqs = new ArrayList<>();
        for (int i = 0; i < TAXA.length; i++) {
            Sequence s = new Sequence();
            s.initByName("taxon", TAXA[i], "value", rows[i]);
            seqs.add(s);
        }
        StandardData dt = new StandardData();
        dt.initByName("nrOfStates", nStates);
        PerCharacterAscertainedAlignment aln = new PerCharacterAscertainedAlignment();
        aln.initByName("sequence", seqs, "dataType", "standard", "userDataType", dt,
                "ascertained", true, "exchangeableStates", exchangeable);
        return aln;
    }

    double logP(TreeLikelihood lik, Alignment aln, StandardData dt, boolean ordered) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", aln, "newick", NEWICK, "IsLabelledNewick", true, "adjustTipHeights", false);
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", substModel(dt, ordered));
        lik.initByName("data", aln, "tree", tree, "siteModel", siteModel);
        return lik.calculateLogP();
    }

    SubstitutionModel substModel(StandardData dt, boolean ordered) throws Exception {
        if (!ordered) {
            LewisMK mk = new LewisMK();
            mk.initByName("datatype", dt);
            return mk;
        }
        int nStates = dt.getStateCount();
        StringBuilder equal = new StringBuilder();
        for (int i = 0; i < nStates; i++) equal.append(1.0 / nStates).append(" ");
        RealParameter p = new RealParameter();
        p.initByName("value", equal.toString().trim());
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", p);
        Ordinal o = new Ordinal();
        o.initByName("stateNumber", nStates, "frequencies", freqs);
        return o;
    }
}
