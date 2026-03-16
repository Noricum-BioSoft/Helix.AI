# Statistical Guidelines and Interpretation Guardrails

## Overview

This document defines statistical validity requirements and interpretation guardrails for all bioinformatics analyses in Helix.AI. These guidelines ensure scientific rigor, reproducibility, and accurate interpretation of results.

---

## Core Principles

### 1. Transparency
- State all assumptions explicitly
- Report methods and parameters clearly
- Disclose limitations and caveats

### 2. Rigor
- Use appropriate statistical methods
- Apply multiple-testing correction by default
- Report effect sizes alongside p-values

### 3. Reproducibility
- Document all parameters and versions
- Report random seeds where applicable
- Provide complete provenance

---

## Statistical Reporting Requirements

### Effect Sizes and Uncertainty

**Requirement:** Always report effect sizes with confidence intervals or standard errors.

**Why:** P-values alone don't indicate biological significance. Effect sizes show the magnitude of differences.

**Examples:**

✅ **Good:**
```
Gene X is upregulated in treatment vs. control:
- Log2 fold change: 2.5 (95% CI: 1.8 to 3.2)
- P-value: 1.2e-10
- FDR: 3.4e-8
```

❌ **Bad:**
```
Gene X is upregulated (p < 0.001)
```

**For common analyses:**

| Analysis Type | Effect Size | Uncertainty Measure |
|---|---|---|
| Differential expression | Log2 fold change | Standard error or 95% CI |
| GWAS | Odds ratio / Beta | 95% CI |
| Survival analysis | Hazard ratio | 95% CI |
| Clustering | Silhouette score | Standard deviation across runs |
| Enrichment | Fold enrichment | 95% CI or p-value |

---

### Multiple Testing Correction

**Requirement:** Apply multiple-testing correction by default for all tests involving multiple hypotheses.

**Why:** Without correction, false positives accumulate rapidly. If testing 20,000 genes at α=0.05, expect 1,000 false positives by chance alone.

**Default Methods:**

| Analysis Type | Default Correction | Alternative |
|---|---|---|
| Differential expression | Benjamini-Hochberg (FDR) | Bonferroni (if few tests) |
| Enrichment analysis | Benjamini-Hochberg (FDR) | Bonferroni |
| GWAS | Bonferroni | FDR |
| Variant calling | Bonferroni or FDR | -- |

**Always state which correction was used:**

✅ **Good:**
```
"After Benjamini-Hochberg FDR correction (FDR < 0.05), 1,234 genes were differentially expressed."
```

❌ **Bad:**
```
"1,234 genes were differentially expressed (p < 0.05)."
```

---

### Sample Size and Statistical Power

**Requirement:** Warn users when sample sizes are insufficient for reliable inference.

**Minimum sample sizes (rules of thumb):**

| Analysis Type | Minimum per Group | Recommended |
|---|---|---|
| Bulk RNA-seq DE | 3 biological replicates | 5-6 replicates |
| scRNA-seq | 100 cells per condition | 1,000+ cells |
| GWAS | 1,000 samples | 10,000+ samples |
| ChIP-seq | 2 replicates | 3+ replicates |
| Proteomics | 3 replicates | 5+ replicates |

**Example warning:**
```
⚠️ Warning: Only 2 replicates per group detected. Results may be unreliable.
Recommendation: Use 3+ biological replicates for differential expression analysis.
```

---

## Analysis-Specific Guidelines

### Differential Expression (RNA-seq)

**Requirements:**

1. **Normalization:**
   - Must normalize for library size and composition
   - State normalization method (TMM, DESeq2, CPM, etc.)

2. **Filtering:**
   - Filter low-count genes before testing
   - State filtering threshold (e.g., "genes with <10 total counts across all samples")

3. **Model specification:**
   - Account for experimental design (blocking factors, paired samples)
   - Warn about confounding variables

4. **Effect size thresholds:**
   - Report both statistical significance (FDR) and biological significance (fold change)
   - Suggest thresholds: FDR < 0.05 AND |log2FC| > 1

**Example output:**
```json
{
  "summary": "Found 1,234 differentially expressed genes (FDR < 0.05, |log2FC| > 1)",
  "methods": {
    "normalization": "DESeq2 median-of-ratios",
    "filtering": "Genes with <10 total counts across all samples removed",
    "model": "~condition + batch",
    "correction": "Benjamini-Hochberg FDR",
    "thresholds": {
      "fdr": 0.05,
      "log2fc_min": 1.0
    }
  },
  "warnings": [
    "Sample sizes are small (n=3 per group). Consider validation with larger cohort."
  ]
}
```

---

### Enrichment Analysis

**Requirements:**

1. **Background set:**
   - State which genes were used as background
   - Default: all genes expressed in the dataset (not all genes in genome)

2. **Multiple testing:**
   - Apply FDR correction across all tested pathways/GO terms

3. **Effect size:**
   - Report fold enrichment, not just p-value
   - Report number of genes in overlap

**Example output:**
```json
{
  "pathway": "Cell cycle (GO:0007049)",
  "genes_in_pathway": 250,
  "genes_in_de_list": 1234,
  "overlap": 45,
  "expected": 15.4,
  "fold_enrichment": 2.92,
  "p_value": 1.2e-8,
  "fdr": 3.4e-6
}
```

---

### Clustering Analysis

**Requirements:**

1. **Cluster validation:**
   - Report silhouette score or other validation metric
   - Warn if clusters are poorly separated

2. **Reproducibility:**
   - Set random seed for stochastic methods (k-means, etc.)
   - Report seed in provenance

3. **Parameter selection:**
   - Justify number of clusters (elbow plot, silhouette analysis)
   - Try multiple resolutions for hierarchical clustering

**Example output:**
```json
{
  "method": "Leiden clustering",
  "resolution": 1.0,
  "n_clusters": 12,
  "silhouette_score": 0.62,
  "modularity": 0.83,
  "warnings": [
    "Cluster 7 has low silhouette score (0.35), may represent heterogeneous cells"
  ]
}
```

---

### Principal Component Analysis (PCA)

**Requirements:**

1. **Variance explained:**
   - Report % variance explained by each PC
   - Plot scree plot or cumulative variance

2. **Confounders:**
   - Check if PCs correlate with batch, sequencing depth, or other technical factors
   - Warn if batch effects are visible

**Example output:**
```json
{
  "pc1_variance": 52.3,
  "pc2_variance": 21.4,
  "total_variance_top10": 88.7,
  "warnings": [
    "PC1 strongly correlates with sequencing batch (r=0.85). Consider batch correction."
  ]
}
```

---

## Assumptions and Confounders

### Always Surface Key Assumptions

**Examples:**

1. **Differential expression:**
   - "Assumed negative binomial distribution for count data"
   - "Assumed independent biological replicates"
   - "Assumed no systematic batch effects"

2. **Variant calling:**
   - "Assumed diploid genome"
   - "Assumed uniform coverage (no CNVs)"
   - "Used default quality filters (QUAL > 30)"

3. **Phylogenetics:**
   - "Assumed molecular clock (constant evolutionary rate)"
   - "Assumed no recombination"
   - "Used GTR+G substitution model"

### Always Flag Potential Confounders

**Common confounders:**

| Confounder | How to Detect | What to Report |
|---|---|---|
| Batch effects | PCA, clustering | "Samples cluster by batch rather than condition" |
| Sequencing depth | Correlation with PC1 | "PC1 correlates with library size (r=0.75)" |
| Sample quality | QC metrics | "Sample X has low quality (mean Q-score = 28)" |
| Outliers | PCA, distance metrics | "Sample Y is an outlier (>3 SD from mean)" |
| Low sample size | Count replicates | "Only 2 replicates per group, results may be unstable" |

---

## Diagnostic Plots

### Required Diagnostics

**For differential expression:**
- ✅ PCA plot colored by condition
- ✅ PCA plot colored by batch (if applicable)
- ✅ MA plot (log2FC vs. mean expression)
- ✅ Volcano plot (log2FC vs. -log10 p-value)
- ✅ Dispersion plot (if using DESeq2/edgeR)

**For scRNA-seq:**
- ✅ QC metrics violin plots (nGenes, nUMIs, % mitochondrial)
- ✅ UMAP/tSNE colored by cluster
- ✅ UMAP/tSNE colored by QC metrics (to check for confounding)
- ✅ Marker gene heatmap

**For clustering:**
- ✅ Silhouette plot
- ✅ Elbow plot (if using k-means)
- ✅ Dendrogram (if hierarchical)

---

## Reporting Statistical Results

### Standard Format

For each statistical test, report:

1. **Test performed:** Name of test (t-test, DESeq2, etc.)
2. **Null hypothesis:** What is being tested
3. **Test statistic:** Value of the statistic
4. **P-value:** Raw p-value
5. **Adjusted p-value:** After multiple testing correction
6. **Effect size:** Magnitude of effect with CI
7. **Sample size:** n for each group

### Example Report

```markdown
## Differential Expression Results

**Method:** DESeq2 v1.34.0 with Benjamini-Hochberg FDR correction

**Null hypothesis:** No difference in gene expression between treatment and control

**Results:**
- Tested: 15,432 genes (after filtering genes with <10 total counts)
- Significant: 1,234 genes (FDR < 0.05)
- Upregulated: 645 genes (log2FC > 1, FDR < 0.05)
- Downregulated: 589 genes (log2FC < -1, FDR < 0.05)

**Sample sizes:**
- Control: n=3 biological replicates
- Treatment: n=3 biological replicates

**Top result:**
- Gene: GENE1
- Base mean: 1234.5
- Log2 fold change: 3.2 (SE: 0.2)
- P-value: 1.2e-50
- FDR: 2.4e-46

**Warnings:**
- Small sample sizes (n=3 per group) may limit statistical power
- Consider validation with qRT-PCR or larger cohort
```

---

## Common Pitfalls to Avoid

### ❌ Don't

1. **Don't report only p-values without effect sizes**
   - P-values depend on sample size; tiny effects can be "significant" with large n

2. **Don't ignore multiple testing**
   - Without correction, most "significant" results are false positives

3. **Don't confuse statistical and biological significance**
   - FDR < 0.05 but log2FC = 0.1 is statistically significant but biologically meaningless

4. **Don't assume normality without checking**
   - Count data is not normally distributed; use appropriate tests (DESeq2, edgeR, etc.)

5. **Don't ignore batch effects**
   - Check PCA plots colored by batch; correct if necessary

6. **Don't cherry-pick results**
   - Report all tests performed, not just "significant" ones

7. **Don't over-interpret exploratory analyses**
   - PCA, clustering, and heatmaps are exploratory; follow up with hypothesis testing

### ✅ Do

1. **Do report effect sizes with confidence intervals**
2. **Do apply multiple testing correction**
3. **Do state all assumptions explicitly**
4. **Do provide diagnostic plots**
5. **Do check for confounders**
6. **Do report negative results** ("no significant difference found")
7. **Do provide complete provenance** (tools, versions, parameters)

---

## Validation and Reproducibility

### Validation Requirements

**For high-impact findings:**
- ✅ Validate top hits with orthogonal method (qRT-PCR for DE genes)
- ✅ Check in independent cohort if available
- ✅ Verify that results are robust to parameter choices

**For exploratory analyses:**
- ✅ Try multiple parameter settings
- ✅ Check stability across random seeds (for stochastic methods)
- ✅ Report sensitivity to outlier removal

### Reproducibility Checklist

- [ ] All random seeds documented
- [ ] All software versions recorded
- [ ] All parameters explicitly stated
- [ ] All data processing steps logged
- [ ] All filtering criteria documented
- [ ] All assumptions stated
- [ ] All confounders checked

---

## Related Documents

- `docs/OUTPUT_SCHEMA.md` - Output format for statistical results
- `docs/TASK_ARCHITECTURE.md` - Workflow composition and provenance
- `agents/workflow-planner-agent.md` - Statistical analysis workflows

---

## References

1. **Multiple testing correction:**
   - Benjamini & Hochberg (1995). Controlling the false discovery rate. J Royal Stat Soc B.

2. **Differential expression:**
   - Love, Huber & Anders (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology.

3. **Effect sizes:**
   - Sullivan & Feinn (2012). Using effect size—or why the P value is not enough. J Grad Med Educ.

4. **Statistical power:**
   - Schurch et al. (2016). How many biological replicates are needed in an RNA-seq experiment. RNA.
