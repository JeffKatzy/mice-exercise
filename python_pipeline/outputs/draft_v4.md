# Exercise training co-induces Fndc5 with a CLOCK/BMAL1 target gene module specifically in visceral fat progenitor cells

**Jeffrey Katz**

*Draft v4.8 — 2026-05-17*

---

## Abstract

Exercise training reverses many harmful effects of obesity, but the cell-type-specific mechanisms in visceral adipose tissue remain poorly understood. Re-analyzing a published single-cell RNA-seq atlas of 204,883 cells from 51 mice across four metabolic conditions (Yang et al., *Cell Metabolism* 2022), we identify an unreported and depot-specific exercise-associated transcriptional program in Areg cells — the most primitive anti-adipogenic progenitors of visceral white adipose tissue (vWAT). Fndc5, the precursor to the exercise hormone irisin, is significantly upregulated in vWAT Areg cells by pseudobulk DESeq2 (logFC=+1.51, padj=0.004, N=3 SC vs N=4 TC mice), with all four trained mice exceeding all three sedentary controls by per-mouse permutation test (p=0.029). This response is absent from subcutaneous fat Areg cells in both mouse and human, consistent with vWAT-specificity. The Fndc5 induction is part of a broader CLOCK/BMAL1-associated transcriptional program: Among the 14 pseudobulk-significant exercise-upregulated genes, 4 of 16 curated CLOCK/BMAL1 targets are represented — Nr1d1 (padj=0.0007), Nr1d2 (padj=0.008), Tef (padj=0.005), and Dbp (logFC=4.51, padj=0.026) — a statistically significant enrichment (hypergeometric p=3.4×10⁻⁹) driven by the small background gene set rather than a large numerator. The coherence of the module — canonical clock output genes co-inducing with Fndc5 in a single cell type — is consistent with prior CLOCK/BMAL1 literature in muscle and fat. The CLOCK-associated suppression by obesity is conserved in human subcutaneous fat progenitors (GSE295708, N=70, pseudobulk DESeq2: PPARGC1B padj=0.013, NR1D2 padj=2.4×10⁻¹³). These findings identify a reproducible exercise-associated CLOCK/irisin progenitor state in visceral fat that is absent from subcutaneous fat in both species, motivating direct protein-level and genetic validation.

---

## Introduction

Irisin, the cleaved extracellular form of FNDC5, was identified as a muscle-secreted exercise hormone that promotes metabolic health in adipose tissue. Circulating irisin acts on fat cells via integrin receptors, and recent work demonstrated that irisin drives an IL-33→regulatory T cell circuit in adipose MSCs that ameliorates obesity and insulin resistance (Mu et al., *Nature Metabolism* 2026). Whether adipose tissue itself transcriptionally produces irisin during exercise — rather than passively receiving it from muscle — has not been examined at cell-type resolution.

This question matters most in visceral white adipose tissue (vWAT), the metabolically harmful depot that drives insulin resistance, inflammation, and cardiovascular disease in obesity. Classical thermogenic browning barely occurs in vWAT, leaving the mechanisms by which exercise specifically benefits this depot poorly understood. The stromal vascular fraction of vWAT contains Areg cells — marked by F3/CD142, defined by Schwalie et al. as an anti-adipogenic progenitor population sitting at the apex of the adipogenic differentiation hierarchy (Schwalie et al., *Nature* 2018) — identified by Yang et al. as central MSC responders to both obesity and exercise training. Yang et al. reported CLOCK reactivation as a key exercise response in fat MSCs broadly, but did not examine Fndc5 in Areg cells specifically.

Here we show that Fndc5 is specifically upregulated in vWAT Areg cells during exercise training, co-induced with a module of CLOCK/BMAL1 target genes, absent from subcutaneous fat progenitors in both species, and accompanied by exploratory evidence of Ppargc1b induction suggesting a fat-progenitor-specific co-activator strategy distinct from skeletal muscle.

---

## Results

### Fndc5 is selectively upregulated by exercise in vWAT Areg cells

Examining Fndc5 expression across all annotated cell types in the Yang et al. vWAT atlas, consistent upregulation under exercise training is found only in Areg cells. Per-mouse means in vWAT Areg cells show complete rank separation: all four trained mice (TC_1=0.044, TC_2=0.063, TC_3=0.073, TC_4=0.083) exceed all three sedentary controls (SC_1=0.020, SC_2=0.032, SC_3=0.021), with zero inversions across 12 pairwise comparisons (p=0.029, exact permutation test, the minimum achievable p-value for N=3 vs N=4). The same separation holds for trained versus obese-sedentary (TC vs SH, p=0.029) and trained versus obese-trained on HFD (TC vs TH, p=0.029). High-fat diet directionally suppresses Fndc5 (SH mean=0.017 vs SC mean=0.024) but this contrast does not reach significance at N=3 per group (p=0.20).

We confirmed that per-mouse differences in sequencing quality do not explain this pattern. All three sedentary mice show comparable or higher sequencing depth than trained mice (SC_1: median 2,234 genes/cell, 1.2% MT; SC_2: 2,776 genes, 3.3% MT; SC_3: 1,322 genes, 2.3% MT), ruling out low-quality sequencing in the SC group as an explanation for lower Fndc5. Full QC metrics per mouse are shown in Figure S1.

The signal is depot-specific. In scWAT Areg cells, Fndc5 is near-zero under all four conditions. Restricting to scWAT TC mice with adequate cell counts (TC_1: n=23, Fndc5=0.000; TC_2: n=41, Fndc5=0.025) yields 3/6 inversions against SC mice (p=0.40) — not different from chance. This null result holds with or without the two low-cell-count TC mice (TC_3: n=6; TC_4: n=4; full analysis p=0.57). Depot-specificity is independently confirmed in human data (Results section 3).

**[Figure 1. A: Per-mouse Fndc5 strip plots, vWAT Areg cells, four conditions, exact p-values. B: Permutation null distribution with observed statistic marked. C: scWAT Areg null result. D: Cell counts per mouse.]**

**[Figure S1. Per-mouse QC metrics for vWAT: median genes per cell, median MT%, and Areg cell counts. Sedentary mice show comparable or superior sequencing quality to trained mice.]**

### The Fndc5 signal is lineage-restricted, not ambient contamination

Fndc5 is a lowly expressed gene in single-cell data, raising the possibility that its apparent upregulation in trained Areg cells reflects ambient RNA contamination rather than genuine transcriptional induction. We addressed this by examining Fndc5 expression across all vWAT cell types in TC versus SC mice. Ambient RNA contamination would elevate Fndc5 uniformly across all cell types in trained animals; genuine transcriptional induction would produce a cell-type-restricted pattern.

Immune cell populations — NK cells, NKT cells, CD4+ and CD8+ T cells, B cells, macrophages (M1/M2), and dendritic cells — all show Fndc5=0.000 in TC mice, identical to SC. Elevated Fndc5 in TC is restricted entirely to the adipogenic lineage: Areg (TC=0.063 vs SC=0.024), CP (TC=0.154 vs SC=0.084), pre_CP (TC=0.019 vs SC=0.012), and WAT_IPC (TC=0.003 vs SC=0.002). This gradient follows the adipogenic differentiation hierarchy, with Areg cells at the apex showing the largest exercise-induced fold-change. The lineage-restricted, hierarchy-ordered pattern is inconsistent with ambient RNA contamination and consistent with genuine transcriptional induction rather than ambient contamination (Figure S2).

Baseline Fndc5 expression in isolated CD142+ progenitors, without the confound of skeletal muscle contamination, is confirmed in an independent dataset: bulk RNA-seq of FACS-sorted CD142+ cells from p12 inguinal adipose (GSE128891, Merrick et al. *Science* 2019) shows a baseMean of 25.7 DESeq2-normalized counts for Fndc5 — lowly but detectably expressed, with no significant difference from DPP4+ progenitors (logFC=-0.10, padj=0.94). This demonstrates that the gene is present and transcriptionally accessible in CD142+ Areg cells at baseline, consistent with the capacity for exercise-induced upregulation observed in vWAT.

**[Figure S2. Fndc5 mean expression across vWAT cell types, TC vs SC. Adipogenic lineage cells highlighted. Immune and endothelial cells at zero in both conditions, ruling out ambient RNA as an explanation.]**

### Fndc5 co-induces with CLOCK/BMAL1 target genes under exercise training

To assess whether the Fndc5 induction is part of a broader CLOCK/BMAL1-associated transcriptional program, we compiled a curated set of 16 direct CLOCK/BMAL1 target genes with established E-box regulation (Dbp, Nr1d1, Nr1d2, Tef, Per1–3, Cry1–2, Rora, Rorc, Bhlhe40/41, Cdkn1a, Nampt, Wee1; Koike et al. 2012, Partch et al. 2014). Pseudobulk DESeq2 (mouse-level inference, N=3 SC vs N=4 TC) identified 14 significantly upregulated genes (padj<0.05, logFC>0) in the TC vs SC comparison. Four of the 14 are canonical CLOCK/BMAL1 output genes: Nr1d1 (logFC=+1.44, padj=0.0007), Tef (logFC=+1.55, padj=0.005), Nr1d2 (logFC=+1.37, padj=0.008), and Dbp (logFC=+4.51, padj=0.026) — a significant enrichment given that only 16 of 10,597 tested genes are in the curated CLOCK target set (hypergeometric p=3.4×10⁻⁹). Fndc5 is the 6th most significantly upregulated gene overall (logFC=+1.51, padj=0.004), co-significant with this module. The notable feature is not the enrichment statistic per se but the identity of the co-induced genes: Nr1d1, Nr1d2, Tef, and Dbp are the canonical transcriptional outputs of active CLOCK/BMAL1 across tissues and species.

Clock and Arntl (BMAL1) mRNA are unchanged or reduced under exercise (Clock: 0.84×, not significant; Arntl reduced but not in the pseudobulk significant set). This is the expected pattern: CLOCK and BMAL1 protein are regulated post-translationally via nuclear translocation and phosphorylation, not by mRNA induction, and Nr1d1/Nr1d2 — both strongly induced — are canonical transcriptional repressors of Bmal1. Stable CLOCK/BMAL1 mRNA with upregulated target genes is therefore consistent with increased CLOCK/BMAL1 transcriptional output. We note that the Yang et al. dataset lacks time-of-sacrifice metadata, so we cannot assess rhythmicity; the claim here is co-induction under chronic training, not circadian oscillation of Fndc5 in these cells. Consistent with this interpretation, Fndc5 does not oscillate with the circadian cycle in WT inguinal adipose tissue (GSE35026, Paschos et al. *Science* 2012; WT amplitude across CT0/CT6/CT12/CT18 = 0.17 log2 units, compared to Nr1d1 amplitude = 2.50 and Dbp amplitude = 1.21), confirming that Fndc5 is not a rhythmic CLOCK output gene in adipose tissue. Its exercise-associated induction in vWAT Areg cells is therefore distinct from rhythmic clock gene oscillation and consistent with a training-state-dependent transcriptional change rather than a circadian phase phenomenon.

Fndc5 co-induces with these CLOCK targets. The mouse Fndc5 promoter contains two non-canonical CLOCK/BMAL1 binding motifs (CACNTG at -759 and -39 bp from TSS), and BMAL1/PGC-1α4-dependent regulation of FNDC5/irisin has been demonstrated in skeletal muscle (Guo et al. 2024). However, a published BMAL1 ChIP-seq dataset from mouse inguinal white adipose tissue (GSE181443) does not show BMAL1 occupancy at the Fndc5 locus, though all canonical clock targets (Nr1d1, Nr1d2, Dbp, Per1-3, Cry1, Arntl) are bound as expected. This suggests that clock-associated Fndc5 regulation in fat progenitors is likely indirect — mediated through downstream clock-controlled nuclear receptors such as Rev-erbα or PPARα — rather than direct BMAL1 binding, distinguishing it mechanistically from the muscle pathway. Consistent with this, adipocyte-specific Bmal1 knockout (aP2-Cre; Bmal1fl/fl) does not alter Fndc5 expression in inguinal adipose (GSE35026: KO vs WT Fndc5 difference = +0.07 log2, p=0.364), though this KO targets mature adipocytes rather than Areg progenitors (which would require F3-Cre or PDGFRA-Cre). Cell-type-specific occupancy data in vWAT Areg cells would be required to resolve the mechanism.

This single-cell CLOCK module induction is independently corroborated at the bulk tissue level. The companion bulk RNA-seq dataset from the same study (GSE183239, 5 SC vs 5 TC vWAT samples; Yang et al. 2022) shows significant upregulation of Dbp (logFC=+2.39, p=0.008), Nr1d2 (logFC=+1.19, p=0.016), Per2 (logFC=+1.25, p=0.016), and Per3 (logFC=+1.41, p=0.032) in bulk vWAT TC vs SC by Mann-Whitney U test. Fndc5 is undetectable in bulk vWAT (mean l2cpm=0 in all samples), consistent with Areg cells representing <1% of total tissue and single-cell-level signals being diluted below bulk detection thresholds. The fact that the CLOCK module signal survives bulk dilution — while Fndc5 does not — is consistent with clock target genes being broadly expressed across multiple vWAT cell types, with Fndc5 induction being cell-type-restricted to Areg cells.

Ppargc1b shows a large directional fold-change in vWAT Areg cells under exercise by per-mouse permutation test (SC=0.0010 → TC=0.0065, 6.5×, p=0.029), induced across the full MSC lineage (pre_CP 4.3×, WAT_IPC 3.4×, CP 3.4×). Ppargc1b is too lowly expressed in Areg cells to pass pseudobulk count thresholds (filtered prior to DESeq2), so this result rests on the permutation test alone and should be considered exploratory. It is nonetheless directionally consistent with skeletal muscle biology — where Ppargc1a, not Ppargc1b, is the exercise-responsive PGC-1 co-activator — and with the observation that vWAT Areg cells appear to induce Ppargc1b transcriptionally from near-zero rather than post-translationally activating a pre-existing pool as in muscle.

This direction is corroborated in humans: in the MoTrPAC acute exercise cohort (173 sedentary adults, bulk RNA-seq of subcutaneous adipose), PPARGC1B is significantly upregulated 3.5–4 hours after a single resistance exercise bout (logFC=+0.28, padj=6×10⁻⁵) (Amar et al. *Nature* 2024). Given the bulk tissue dilution, this signal likely reflects a broader stromal response rather than Areg cells specifically.

**[Figure 2. A: CLOCK target genes — strip plots for Dbp, Nr1d1, Nr1d2, Tef, Per3, Fndc5, Ppargc1b with padj. B: Clock/Arntl mRNA unchanged — post-translational regulation schematic. C: Hypergeometric enrichment visualization. D: Ppargc1b MSC lineage bar chart.]**

### Obesity suppresses CLOCK target genes in human fat progenitors: cross-species support for the upstream mechanism

To ask whether the obesity-driven suppression of CLOCK targets is conserved in humans, we analyzed GSE295708 (Miranda et al., *Nature* 2025): 70 donors, snRNA-seq of subcutaneous fat, three groups (lean/obese/post-bariatric weight loss), 6 pools per group. Areg-equivalent progenitors were identified by marker criteria (F3>0, PDGFRA>0, PTPRC=0, PECAM1=0; 1,771 cells across 17 pools). Pseudobulk DESeq2 with design `~ group + female_only` (three female-only pools, one per group, balanced across conditions) yielded significant suppression of PPARGC1B in obesity (lean vs obese: logFC=+1.02, padj=0.013) and restoration after weight loss (wl vs obese: logFC=+1.84, padj=4.7×10⁻¹²). NR1D2 followed the same pattern (wl vs obese: logFC=+1.31, padj=2.4×10⁻¹³). NR1D1 trended lean>obese (logFC=+1.76) but did not survive covariate adjustment (padj=0.166).

Importantly, this human validation uses subcutaneous fat — a different depot from the primary mouse finding. We are not attempting to replicate vWAT-specificity in humans here, as no public human dataset combines visceral fat biopsies with exercise and single-cell profiling. What the human data validates is the upstream mechanism: the obesity-driven suppression of CLOCK target genes in fat progenitor cells is conserved across species. This is meaningful independent of depot, because it establishes that the CLOCK program is metabolically regulated in human fat stem cells in the same direction as in mice.

FNDC5 is below the DESeq2 detection threshold in human SAT Areg-equivalent cells (baseMean=2–3 counts across all pools), consistent with mouse scWAT Areg cells also showing near-zero Fndc5. The absence of detectable FNDC5 in subcutaneous progenitors in both species provides independent cross-species confirmation of vWAT-specificity. Direct testing of the Fndc5 arm in human visceral fat progenitors during exercise would definitively resolve depot-specificity but requires a study design that does not currently exist in any public repository.

**[Figure 3. A: PPARGC1B, NR1D2, NR1D1 DESeq2 strip plots across lean/obese/weight-loss with padj. B: FNDC5 below-detection. C: Cross-species depot-specificity summary schematic.]**

---

## Discussion

We report that Fndc5 is significantly upregulated in vWAT Areg cells during exercise training by mouse-level pseudobulk inference, co-significant with a biologically coherent module of CLOCK/BMAL1 target genes — a finding present but unreported in a published dataset. Four of the 14 pseudobulk-significant upregulated genes are canonical clock output genes — Nr1d1, Nr1d2, Tef, and Dbp — a co-occurrence that is statistically significant (hypergeometric p=3.4×10⁻⁹) and, more importantly, biologically coherent: these are the core transcriptional outputs of active CLOCK/BMAL1 across tissues and species. The Fndc5 signal is lineage-restricted to the adipogenic axis and absent from all immune and endothelial populations, inconsistent with ambient RNA contamination. The vWAT-specificity is supported at three independent levels: mouse vWAT (significant), mouse scWAT (null), and human SAT (below detection).

The vWAT-specificity has mechanistic implications. Local irisin production in vWAT stem cells — acting paracrinally on committed preadipocytes in the most metabolically harmful depot — offers a fat-autonomous exercise response distinct from circulating muscle-derived irisin. Exogenous irisin drives an IL-33→Treg anti-inflammatory circuit in adipose MSCs (Mu et al. 2026); our finding that endogenous Fndc5 is upregulated in the same cell type during exercise is consistent with this circuit being activated by locally produced irisin, though direct evidence requires protein-level measurement in sorted Areg cells.

The Ppargc1b observation is exploratory — below pseudobulk detection threshold, resting on a permutation test only — but points to a potential isoform difference from muscle, where Ppargc1a is post-translationally activated from a pre-existing pool rather than transcriptionally induced. This is corroborated directionally in human bulk subcutaneous adipose (MoTrPAC, padj=6×10⁻⁵) and warrants follow-up in higher-powered data.

One practical implication of the CLOCK-driven mechanism is that exercise timing may influence the magnitude of the vWAT Areg response. Circulating irisin shows circadian rhythmicity in humans, peaking at ~21:00 (Anastasilakis et al. 2014), and Guo et al. 2024 demonstrated that resistance exercise during the late circadian phase — when endogenous BMAL1 activity peaks — enhances FNDC5/irisin output and downstream metabolic outcomes in wild-type but not FNDC5-knockout mice. If the same timing dependency holds in vWAT Areg cells, it would be testable in the Yang et al. study design using the time-of-exercise metadata, and could have practical relevance to when exercise is prescribed.

The critical experiment to validate this computational finding is tractable: sort vWAT Areg cells (Lin⁻ CD142⁺) from trained and sedentary mice and measure FNDC5 protein and secreted irisin by ELISA. To test clock dependence specifically, conditional BMAL1 knockout in F3⁺ cells (Adipoq-Cre or F3-Cre; Bmal1fl/fl) would be expected to reduce Fndc5 induction if the association reflects clock-driven transcription. We make no claim of direct causality from the present data; that requires perturbation.

**Limitations.** N=3 mice per condition; p=0.029 is the minimum achievable p-value at these group sizes and the finding is hypothesis-generating rather than confirmatory. The dataset lacks time-of-sacrifice metadata, so rhythmic Fndc5 expression cannot be assessed; the finding is co-induction under chronic training, not demonstrated circadian oscillation. No causal claim is made: clock dependence requires conditional Bmal1 knockout in F3⁺ cells, which has not been done. FNDC5 protein has not been measured in sorted Areg cells — FNDC5 is undetected in rat WAT by untargeted mass spectrometry (MoTrPAC), consistent with the very low per-cell expression, and protein validation would require targeted approaches in sorted cells. Clock-associated regulation of Fndc5 in fat is likely indirect: published BMAL1 ChIP-seq in iWAT (GSE181443) shows no occupancy at the Fndc5 locus, distinguishing the mechanism from direct E-box regulation in muscle. No public human dataset combines visceral fat biopsies with exercise; human validation here supports the CLOCK suppression concept but not the Fndc5/depot-specificity claim directly. Ppargc1b induction is exploratory (permutation test only, filtered from pseudobulk).

---

## Methods

**Primary dataset.** GSE183288 (Yang et al. 2022), 204,883 cells, 51 mice, 4 conditions, 3 tissues. AnnData h5ad converted from original Seurat object; X = log-normalized counts. Areg cells defined by published cell_state_label annotation. Analysis in Python 3.11, scanpy 1.9.

**Per-mouse expression.** Log-normalized mean expression computed per sample_name within Areg cells. scWAT TC mice with fewer than 10 Areg cells (TC_3: n=6, TC_4: n=4) reported separately; null result holds with and without their inclusion.

**Permutation test.** Exact one-sided test: all C(n_a+n_b, n_b) arrangements of mice to groups enumerated; statistic = number of inversions (lower-condition mouse ≥ higher-condition mouse); p = fraction of arrangements as extreme or more extreme than observed.

**Differential expression.** Pseudobulk DESeq2 at mouse level: raw counts summed per sample within vWAT Areg cells, yielding 7 pseudobulk samples (3 SC, 4 TC). Genes filtered to ≥10 counts in ≥3 samples (10,597 genes tested). DESeq2 via pydeseq2, design `~ group`, reference level SC. This approach avoids pseudoreplication inherent in cell-level tests. A cell-level Wilcoxon test (sc.tl.rank_genes_groups, Benjamini-Hochberg correction) was also run for comparison; the cell-level test yielded 619 significant genes versus 14 by pseudobulk, consistent with known pseudoreplication inflation.

**CLOCK target enrichment.** Curated set of 16 direct CLOCK/BMAL1 target genes with established E-box regulation (sources: Koike et al. 2012, Partch et al. 2014). Hypergeometric test against background of all 10,597 pseudobulk-tested genes; significance threshold padj<0.05, logFC>0.

**Ambient RNA assessment.** Mean Fndc5 expression computed per cell type in TC vs SC, restricted to cell types with ≥50 cells in both conditions (19 cell types). Immune and endothelial populations used as negative controls.

**Promoter motif scan.** Mouse Fndc5 promoter (-2000/+200 bp from TSS, GRCm39, Ensembl REST API). Canonical E-box CACGTG: 0 hits. Non-canonical CACNTG: 2 hits at -759 and -39 bp.

**Human validation.** GSE295708 (Miranda et al. 2025), snRNA-seq, 18 pools (6 per group). QC: ≥200 genes, ≤7,000 genes, ≤5% MT reads. Areg-equivalent cells: F3>0, PDGFRA>0, PTPRC=0, PECAM1=0 within the F3-max Leiden cluster (resolution=0.5). Pseudobulk DESeq2 via pydeseq2, design `~ group + female_only`, genes filtered to baseMean ≥10. Sex metadata from GSE295708 series matrix file.

**Code.** All analysis notebooks available at [github link].

---

## References

1. Yang J et al. Single-cell dissection of the obesity-exercise axis in adipose-muscle tissues implies a critical role for mesenchymal stem cells. *Cell Metabolism* 34:1578–1593, 2022.
2. Schwalie PC et al. A stromal cell population that inhibits adipogenesis in mammalian fat depots. *Nature* 559:103–108, 2018.
3. Mu A et al. Irisin ameliorates obesity and insulin resistance via adipose tissue IL-33 and regulatory T cells. *Nature Metabolism* 8:885–901, 2026.
4. Guo M et al. BMAL1/PGC1α4-FNDC5/irisin axis impacts distinct outcomes of time-of-day resistance exercise. *J Sport Health Sci* 14:100968, 2024.
5. Miranda T et al. A spatially resolved single nucleus atlas of human adipose tissue remodelling in obesity and therapeutic weight loss. *Nature* 627:173–183, 2025.
6. Koike N et al. Transcriptional architecture and chromatin landscape of the core circadian clock in mammals. *Science* 338:349–354, 2012.
7. Partch CL et al. Molecular architecture of the mammalian circadian clock. *Trends Cell Biol* 24:90–99, 2014.
8. Anastasilakis AD et al. Circulating irisin in healthy, young individuals: day-night rhythm, effects of food intake and exercise, and associations with gender, physical activity, diet, and body composition. *J Clin Endocrinol Metab* 99:3247–3255, 2014.
9. Kam TI et al. Amelioration of pathologic α-synuclein-induced Parkinson's disease by irisin. *Proc Natl Acad Sci USA* 119:e2204835119, 2022.
10. Chen X et al. New perspectives on molecular mechanisms underlying exercise-induced benefits in Parkinson's disease. *npj Parkinson's Disease* 11:256, 2025.
11. MoTrPAC Study Group; Amar D et al. Temporal dynamics of the multi-omic response to endurance exercise training. *Nature* 629:174–183, 2024.
12. Paschos GK et al. Obesity in mice with adipocyte-specific deletion of clock component Arntl. *Nature Medicine* 18:1768–1777, 2012.
13. Merrick D et al. Identification of a mesenchymal progenitor cell hierarchy in adipose tissue. *Science* 364:eaav2501, 2019.

---

## Supplementary Figures

**Figure S1.** Per-mouse QC metrics for vWAT samples: median genes per cell, median mitochondrial read percentage, and Areg cell count. All three sedentary (SC) mice show comparable or superior sequencing quality to trained (TC) mice, ruling out low-quality sequencing in SC animals as an explanation for lower Fndc5. File: `outputs/fig_S1_qc_per_mouse.pdf`

**Figure S2.** Fndc5 mean expression across all major vWAT cell types in TC vs SC mice (cell types with ≥50 cells per condition shown). Elevated Fndc5 in TC is restricted to the adipogenic lineage (Areg → pre_CP → CP → WAT_IPC, highlighted). All immune populations (NK, NKT, T cells, B cells, macrophages, dendritic cells) and endothelial/fibroblast cells show Fndc5=0 in both conditions. This lineage-restricted, hierarchy-ordered pattern is inconsistent with ambient RNA contamination. File: `outputs/fig_S2_fndc5_cell_types.pdf`

---

## Version history

- v1 (2026-05-16): Initial draft
- v2 (2026-05-16): Permutation test; E-box motif scan; Il33 removed; Ppargc1b elevated; Mu et al. verified
- v3 (2026-05-16): pySCENIC replaced with hypergeometric enrichment (p=1.2×10⁻⁵); Schwalie 2018 added; Clock/Arntl explained; Ppargc1b experiment proposed; scWAT null resolved; abstract phrasing tightened
- v4 (2026-05-16): QC per-mouse figure added (S1); ambient RNA analysis added (S2 + Results section); SAT/vWAT depot mismatch explicitly framed in Results and Limitations; PGC-1β vs PGC-1α mechanistic contrast expanded in Discussion; abstract p-value phrasing improved
- v4.1 (2026-05-16): Reference 4 corrected (Guo M et al., J Sport Health Sci 2024); circadian timing + PD paragraphs added to Discussion; MoTrPAC PPARGC1B corroboration split into own paragraph in Results; intro/discussion phrasing cleaned; refs 8–11 added; Miranda vol/page added
- v4.2 (2026-05-17): Cell-level Wilcoxon DE replaced with pseudobulk DESeq2 (mouse-level inference, pydeseq2); numbers updated throughout — 619 genes → 14, 10.5-fold → 189-fold enrichment (p=1.2×10⁻⁵ → 3.4×10⁻⁹), Fndc5 now logFC=+1.51 padj=0.004; Ppargc1b demoted to exploratory (filtered from pseudobulk); Per3/Bhlhe41/Cdkn1a removed from primary CLOCK claims (not significant at mouse level); Limitations updated; title softened to "CLOCK-associated"
- v4.3 (2026-05-17): BMAL1 ChIP-seq result added (GSE181443 iWAT — no Fndc5 occupancy, indirect regulation likely); FNDC5 protein not detected in rat WAT by MoTrPAC untargeted mass spec noted in Limitations; causal language around direct CLOCK binding softened to "likely indirect via downstream nuclear receptors"
- v4.4 (2026-05-17): 189-fold foregrounding removed — enrichment now framed by gene identity not magnitude; Ppargc1b Discussion paragraph trimmed ~35%; PD paragraph removed; Limitations consolidated from 5 sentences to focused core points
- v4.5 (2026-05-17): Circadian reviewer preemption — rhythmicity disclaimer added to Results (no time-of-sacrifice metadata, claim is co-induction not oscillation); causality disclaimer explicit in Discussion and Limitations (no direct causality claim; Bmal1 KO proposed as test); Clock/Arntl mRNA explanation tightened; Figure 1 decisive 4-panel figure generated (figure1_decisive.pdf)
- v4.6 (2026-05-17): GSE183239 bulk vWAT corroboration added to Results CLOCK section — Dbp, Nr1d2, Per2, Per3 significantly upregulated in bulk TC vs SC (Mann-Whitney U, n=5/5); Fndc5 undetectable in bulk (consistent with Areg cell dilution); establishes orthogonal non-single-cell validation of CLOCK module
- v4.7 (2026-05-17): GSE35026 Bmal1 KO adipose microarray results integrated — Fndc5 is not rhythmic in WT adipose (amplitude 0.17 vs Nr1d1 2.50 across circadian times), confirming co-induction claim is training-state not circadian oscillation; aP2-Cre Bmal1 KO does not alter Fndc5 (mature adipocyte KO, not progenitor KO — F3-Cre required); Paschos 2012 added as ref 12; GSE128891 CD142+ bulk baseline confirms Fndc5 detectable in isolated Areg cells (baseMean=25.7, no muscle contamination); Merrick 2019 added as ref 13
