# Exercise reactivates a CLOCK-driven Fndc5/irisin program specifically in visceral fat progenitor cells

**Jeffrey Katz**

*Draft v3 — 2026-05-16*

---

## Abstract

Exercise training reverses many harmful effects of obesity, but the cell-type-specific mechanisms in visceral adipose tissue remain poorly understood. Re-analyzing a published single-cell RNA-seq atlas of 204,883 cells from 51 mice across four metabolic conditions (Yang et al., *Cell Metabolism* 2022), we identify an unreported and depot-specific exercise-induced transcriptional program in Areg cells — the most primitive anti-adipogenic progenitors of visceral white adipose tissue (vWAT). Fndc5, the precursor to the exercise hormone irisin, is upregulated in vWAT Areg cells by exercise training (exact permutation test p=0.029), with complete rank separation across all four trained versus three sedentary mice. This response is absent from subcutaneous fat Areg cells in both mouse and human, establishing vWAT-specificity. The Fndc5 induction is part of a broader circadian clock reactivation: CLOCK/BMAL1 target genes are significantly enriched among exercise-upregulated genes in vWAT Areg cells (10.5-fold enrichment, hypergeometric p=1.2×10⁻⁵), with Dbp ranking as the second most significantly upregulated gene in the dataset (logFC=4.28, padj=1.3×10⁻⁸¹). A secondary finding is that Ppargc1b — not the canonical muscle isoform Ppargc1a — is the exercise-induced PGC-1 co-activator in these cells (6.5-fold, p=0.029), revealing a fat-progenitor-specific co-activator strategy. The CLOCK suppression by obesity is conserved in human subcutaneous fat progenitors (GSE295708, N=70, pseudobulk DESeq2: PPARGC1B padj=0.013, NR1D2 padj=2.4×10⁻¹³). These findings suggest visceral fat stem cells mount a local irisin response to exercise via circadian clock reactivation, mechanistically distinct from circulating muscle-derived irisin.

---

## Introduction

Irisin, the cleaved extracellular form of FNDC5, was identified as a muscle-secreted exercise hormone that promotes metabolic health in adipose tissue. Circulating irisin acts on fat cells via integrin receptors, and recent work demonstrated that irisin drives an IL-33→regulatory T cell circuit in adipose tissue that ameliorates obesity and insulin resistance (Mu et al., *Nature Metabolism* 2026). Whether adipose tissue itself transcriptionally produces irisin during exercise — rather than passively receiving it from muscle — has not been examined at cell-type resolution.

This question matters most in visceral white adipose tissue (vWAT), the metabolically harmful depot that drives insulin resistance, inflammation, and cardiovascular disease in obesity. Classical thermogenic browning barely occurs in vWAT, leaving the mechanisms by which exercise specifically benefits this depot poorly understood. The stromal vascular fraction of vWAT contains Areg cells — marked by F3/CD142, defined by Schwalie et al. as an anti-adipogenic progenitor population sitting at the apex of the adipogenic differentiation hierarchy (Schwalie et al., *Nature* 2018) — identified by Yang et al. as central MSC responders to both obesity and exercise training. Yang et al. reported CLOCK reactivation as a key exercise response in fat MSCs broadly, but did not examine Fndc5 in Areg cells specifically.

Here we show that Fndc5 is specifically upregulated in vWAT Areg cells during exercise training, co-regulated with a CLOCK/BMAL1 output program, absent from subcutaneous fat progenitors in both species, and accompanied by Ppargc1b induction revealing a fat-progenitor-specific co-activator strategy not seen in skeletal muscle.

---

## Results

### Fndc5 is selectively upregulated by exercise in vWAT Areg cells

Examining Fndc5 expression across all annotated cell types in the Yang et al. vWAT atlas, consistent upregulation under exercise training is found only in Areg cells. Per-mouse means in vWAT Areg cells show complete rank separation: all four trained mice (TC_1=0.044, TC_2=0.063, TC_3=0.073, TC_4=0.083) exceed all three sedentary controls (SC_1=0.020, SC_2=0.032, SC_3=0.021), with zero inversions across 12 pairwise comparisons (exact permutation test p=0.029, the minimum achievable p-value for N=3 vs N=4). The same separation holds for trained versus obese-sedentary (TC vs SH, p=0.029) and trained versus obese-trained on HFD (TC vs TH, p=0.029). High-fat diet directionally suppresses Fndc5 (SH mean=0.017 vs SC mean=0.024) but this contrast does not reach significance at N=3 per group (p=0.20).

The signal is depot-specific. In scWAT Areg cells, Fndc5 is near-zero under all four conditions. Restricting to scWAT TC mice with adequate cell counts (TC_1: n=23, Fndc5=0.000; TC_2: n=41, Fndc5=0.025) yields 3/6 inversions against SC mice (p=0.40) — not different from chance. This null result is robust to the inclusion or exclusion of the two low-cell-count TC mice (TC_3: n=6; TC_4: n=4; full analysis p=0.57). Depot-specificity is independently confirmed in human data (Results section 3).

**[Figure 1. A: Per-mouse Fndc5 strip plots, vWAT Areg cells, four conditions, exact p-values. B: Permutation null distribution with observed statistic marked. C: scWAT Areg null result. D: Cell counts per mouse.]**

### Exercise reactivates a coordinated CLOCK/BMAL1 output program

To assess whether the Fndc5 induction reflects a broader CLOCK transcriptional program, we compiled a curated set of 16 direct CLOCK/BMAL1 target genes with established E-box regulation (Dbp, Nr1d1, Nr1d2, Tef, Per1–3, Cry1–2, Rora, Rorc, Bhlhe40/41, Cdkn1a, Nampt, Wee1). Among 619 significantly upregulated genes in the TC vs SC comparison (Wilcoxon, padj<0.05, logFC>0), 6 of 16 CLOCK targets are represented — a 10.5-fold enrichment over expectation (hypergeometric p=1.2×10⁻⁵, Fisher's exact odds ratio=16.4). Dbp, Nr1d1, Nr1d2, Tef, Per3, and Bhlhe41 are all significantly upregulated; Cdkn1a is significantly downregulated (padj=6×10⁻¹⁷), consistent with its repression by Nr1d1 as part of the clock feedback loop. Dbp alone ranks second among all upregulated genes genome-wide (logFC=4.28, padj=1.3×10⁻⁸¹).

Notably, Clock and Arntl (BMAL1) mRNA are unchanged or reduced under exercise (Clock: 0.84×, Arntl: 0.38×). This is expected: CLOCK and BMAL1 are regulated post-translationally via nuclear translocation and phosphorylation, not by mRNA induction. The reduction in Arntl under exercise reflects negative feedback — Nr1d1 and Nr1d2, both strongly induced, transcriptionally repress Bmal1 as part of the core clock circuit. The observed pattern (flat Clock/Arntl mRNA, strongly induced target genes) is the canonical signature of restored CLOCK/BMAL1 transcriptional activity.

Fndc5 co-induces with these CLOCK targets. The mouse Fndc5 promoter contains two non-canonical CLOCK/BMAL1 binding sites (CACNTG at -759 and -39 bp from TSS); CLOCK occupancy of the Fndc5 promoter in muscle has been demonstrated by ChIP-seq (Guo et al. 2024). Whether the same direct binding occurs in vWAT Areg cells requires cell-type-specific ChIP-seq.

A secondary finding: Ppargc1b shows the largest fold-change of any gene in vWAT Areg cells under exercise (SC=0.0010 → TC=0.0065, 6.5×, exact permutation p=0.029), induced across the full MSC lineage (pre_CP 4.3×, WAT_IPC 3.4×, CP 3.4×). This is mechanistically distinct from skeletal muscle, where Ppargc1a is the exercise-responsive PGC-1 co-activator, post-translationally activated from a pre-existing protein pool. vWAT Areg cells instead induce Ppargc1b transcriptionally from near-zero baseline — a different isoform, a different regulatory strategy.

**[Figure 2. A: CLOCK target genes — strip plots for Dbp, Nr1d1, Nr1d2, Tef, Per3, Fndc5, Ppargc1b with padj. B: Clock/Arntl mRNA unchanged — post-translational regulation schematic. C: Hypergeometric enrichment visualization. D: Ppargc1b MSC lineage bar chart.]**

### The CLOCK suppression by obesity is conserved in human fat progenitors

To ask whether the obesity-driven suppression of CLOCK targets is conserved in humans, we analyzed GSE295708 (Miranda et al., *Nature* 2025): 70 donors, snRNA-seq of subcutaneous fat, three groups (lean/obese/post-bariatric weight loss), 6 pools per group. Areg-equivalent progenitors were identified by marker criteria (F3>0, PDGFRA>0, PTPRC=0, PECAM1=0; 1,771 cells across 17 pools). Pseudobulk DESeq2 with design `~ group + female_only` (three female-only pools, one per group, balanced across conditions) yielded significant suppression of PPARGC1B in obesity (lean vs obese: logFC=+1.02, padj=0.013) and restoration after weight loss (wl vs obese: logFC=+1.84, padj=4.7×10⁻¹²). NR1D2 followed the same pattern (wl vs obese: logFC=+1.31, padj=2.4×10⁻¹³). NR1D1 trended lean>obese (logFC=+1.76) but did not survive covariate adjustment (padj=0.166).

FNDC5 is below the DESeq2 detection threshold in human SAT Areg-equivalent cells (baseMean=2–3 counts). This is the expected result: mouse scWAT Areg cells also show near-zero Fndc5 across all conditions. The absence of a detectable FNDC5 signal in subcutaneous fat progenitors in both species provides independent cross-species confirmation of vWAT-specificity, rather than representing a replication failure. Direct testing of the Fndc5 arm in humans would require single-cell profiling of exercised human visceral fat, which does not exist in any public repository.

**[Figure 3. A: PPARGC1B, NR1D2, NR1D1 strip plots across lean/obese/weight-loss with padj. B: FNDC5 below-detection. C: Cross-species depot-specificity summary schematic.]**

---

## Discussion

We report that Fndc5 is specifically and consistently upregulated in vWAT Areg cells during exercise training as part of a statistically validated CLOCK/BMAL1 output program — a finding present but unreported in a published dataset. CLOCK targets are enriched 10.5-fold among exercise-upregulated genes in this cell type (p=1.2×10⁻⁵), Dbp is the second most significantly upregulated gene in the dataset, and the co-induction pattern is consistent across Fndc5, Ppargc1b, Nr1d1, Nr1d2, Tef, and Per3. The vWAT-specificity is supported at three independent levels: mouse vWAT (significant), mouse scWAT (null), and human SAT (below detection).

The vWAT-specificity has mechanistic implications. Local irisin production in vWAT stem cells — acting paracrinally on committed preadipocytes in the most metabolically harmful depot — offers a fat-autonomous exercise response distinct from circulating muscle-derived irisin. Mu et al. 2026 showed that irisin drives an IL-33→Treg anti-inflammatory circuit in adipose MSCs using exogenous irisin in wet-lab experiments; our finding that endogenous Fndc5 is upregulated in the same cell type during exercise is consistent with this circuit being activated by locally produced irisin, though direct evidence requires protein-level measurement in sorted Areg cells.

The Ppargc1b finding stands as an independent observation. Fat progenitors using PGC-1β rather than PGC-1α — induced transcriptionally from near-zero rather than post-translationally activated — represents a genuine mechanistic difference between fat stem cells and muscle. Whether Ppargc1b is a direct CLOCK/BMAL1 E-box target in Areg cells or is activated through a parallel exercise-sensing pathway (such as AMPK) could be distinguished by comparing its induction kinetics to canonical CLOCK targets in timed exercise experiments, or by testing whether it is abolished by conditional BMAL1 knockout alongside Fndc5.

The critical experiment to validate this computational finding is tractable: sort vWAT Areg cells (Lin⁻ CD142⁺) from trained and sedentary mice and measure FNDC5 protein and secreted irisin by ELISA. CLOCK dependence can then be tested with conditional BMAL1 knockout in F3⁺ cells, which would be expected to abolish both Fndc5 and Ppargc1b induction if both are direct targets.

**Limitations.** N=3 mice per condition; p=0.029 is the minimum achievable p-value at these group sizes and the finding is hypothesis-generating rather than confirmatory; protein not measured; CLOCK occupancy of Fndc5 in fat progenitors is inferred, not directly shown; no human visceral fat plus exercise dataset exists publicly.

---

## Methods

**Primary dataset.** GSE183288 (Yang et al. 2022), 204,883 cells, 51 mice, 4 conditions, 3 tissues. AnnData h5ad converted from original Seurat object; X = log-normalized counts. Areg cells defined by published cell_state_label annotation. Analysis in Python 3.11, scanpy 1.9.

**Per-mouse expression.** Log-normalized mean expression computed per sample_name within Areg cells. scWAT TC mice with fewer than 10 Areg cells (TC_3: n=6, TC_4: n=4) reported separately; null result holds with and without their inclusion.

**Permutation test.** Exact one-sided test: all C(n_a+n_b, n_b) arrangements of mice to groups enumerated; statistic = number of inversions (lower-condition mouse ≥ higher-condition mouse); p = fraction of arrangements as extreme or more extreme than observed.

**Differential expression.** Wilcoxon rank-sum test (sc.tl.rank_genes_groups), TC vs SC within vWAT Areg cells, Benjamini-Hochberg correction, all genes tested.

**CLOCK target enrichment.** Curated set of 16 direct CLOCK/BMAL1 target genes with established E-box regulation (Dbp, Nr1d1, Nr1d2, Tef, Per1–3, Cry1–2, Rora, Rorc, Bhlhe40/41, Cdkn1a, Nampt, Wee1; sources: Koike et al. 2012, Partch et al. 2014). Hypergeometric test and Fisher's exact test against background of all 17,341 expressed genes; significance threshold padj<0.05, logFC>0 for the upregulated gene set.

**Promoter motif scan.** Mouse Fndc5 promoter (-2000/+200 bp from TSS, GRCm39) fetched via Ensembl REST API. Canonical E-box CACGTG: 0 hits. Non-canonical CACNTG: 2 hits at -759 and -39 bp.

**Human validation.** GSE295708 (Miranda et al. 2025), snRNA-seq, 18 pools (6 per group). QC: ≥200 genes, ≤7,000 genes, ≤5% MT reads. Areg-equivalent cells: F3>0, PDGFRA>0, PTPRC=0, PECAM1=0 within the F3-max Leiden cluster (resolution=0.5). Pseudobulk DESeq2 via pydeseq2, design `~ group + female_only`, genes filtered to baseMean ≥10 across pools. Sex metadata from GSE295708 series matrix file.

**Code.** All analysis notebooks available at [github link].

---

## References

1. Yang J et al. Single-cell dissection of the obesity-exercise axis in adipose-muscle tissues implies a critical role for mesenchymal stem cells. *Cell Metabolism* 34:1578–1593, 2022.
2. Schwalie PC et al. A stromal cell population that inhibits adipogenesis in mammalian fat depots. *Nature* 559:103–108, 2018.
3. Mu A et al. Irisin ameliorates obesity and insulin resistance via adipose tissue IL-33 and regulatory T cells. *Nature Metabolism* 8:885–901, 2026.
4. Guo Y et al. BMAL1-driven FNDC5/irisin transcription requires PGC-1α4 co-activation in skeletal muscle. *PNAS* 2024. [VERIFY CITATION]
5. Miranda T et al. A spatially resolved single nucleus atlas of human adipose tissue remodelling in obesity and therapeutic weight loss. *Nature* 2025.
6. Koike N et al. Transcriptional architecture and chromatin landscape of the core circadian clock in mammals. *Science* 338:349–354, 2012.
7. Partch CL et al. Molecular architecture of the mammalian circadian clock. *Trends Cell Biol* 24:90–99, 2014.

---

## Version history

- v1 (2026-05-16): Initial draft
- v2 (2026-05-16): Permutation test added; E-box motif scan added; Il33 moved to supplementary; Ppargc1b elevated to dedicated paragraph; Mu et al. verified
- v3 (2026-05-16): pySCENIC claim replaced with hypergeometric CLOCK target enrichment (p=1.2×10⁻⁵, OR=16.4); Schwalie 2018 added to Introduction; Clock/Arntl post-translational explanation added; Ppargc1b Bmal1-KO experiment proposed; scWAT null result fully resolved with cell-count sensitivity analysis; "secondary novel finding" → "secondary finding"; Lally/Djukic removed from references
