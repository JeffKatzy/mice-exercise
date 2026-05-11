# Pipeline Summary: Yang et al. Exercise/Obesity Mouse Study (GSE183288)

The pipeline processes single-cell and bulk RNA-seq data from three tissues (scWAT, vWAT, skeletal muscle) across four conditions (standard chow / HFD × sedentary / exercise training) in mice.

---

## scRNA-seq Pipeline

**The central problem:** raw CellRanger output contains three layers of noise — ambient RNA from lysed cells, low-quality barcodes, and doublets — and 42 samples that must be harmonized into a single atlas without erasing biological signal.

```
CellRanger (HPC)
      │
      ▼
Step 1  cellranger_metrics              Aggregate per-sample QC stats → decide which samples to use
      │
      ▼
Step 2  sample_level_diagnostic_plot    Visualize per-cell distributions → set filtering thresholds
      │
      ▼
Step 3  sample_level_processing         Three sequential passes per sample:
        ├── SoupX                       Remove ambient RNA (fixed 20% rate for this study)
        ├── Seurat filter               Drop low-quality cells (genes, UMIs, %mt thresholds)
        └── DoubletFinder               Score and remove doublets (~3.1% expected rate)
      │
      ▼
Step 4  pseudobulk_clustering           Sanity check: do samples cluster by tissue/condition
                                        before integration? Catches swaps and batch surprises.
      │
      ▼
Step 5  sample_integration              Merge all samples → normalize → PCA → FIt-SNE + UMAP
        ├── fancy_tsne                  Best-practice tSNE: PCA init, multi-scale perplexity,
        │                               learning rate = n/12, FIt-SNE for 200k+ cells
        └── Leiden + DBSCAN             Two complementary clustering approaches
      │
      ▼
Step 6  add_metadata                    Join sample-level phenotype table onto every cell
                                        (tissue, diet, exercise, collection day/time)
      │
      ▼
Step 7  analysis_pipeline_scwat_vwat_skm    Project-specific downstream analysis:
        ├── feature_plot                    Render any gene or metadata variable onto tSNE/UMAP
        └── (cell type annotation, DEGs, pathway enrichment, CCC, GRN)
```

### Key Design Decisions

| Decision | Value | Rationale |
|----------|-------|-----------|
| `soupx_mode` | `fixed_0.2` | Conservative ambient correction for fragile adipose/muscle tissue |
| `per_mt_high` (per sample) | `30%` | Adipose and muscle have naturally higher mitochondrial content; permissive first pass |
| `mt_cutoff` (at integration) | `10%` | Stricter second pass catches cells that look atypical in the full atlas context |
| `integration_method` | `merge` | No batch correction by design — diet × exercise signal must not be regressed out; Step 4 confirmed samples cluster by tissue rather than batch |
| `rm_prolif` | `True` | Proliferating cells share a cell-cycle transcriptome across lineages, creating spurious clusters |
| Atlas structure | Full atlas + per-tissue | Full atlas for cross-tissue comparisons; per-tissue objects for higher-resolution within-tissue analysis |

### Study Design

| Code | Diet | Exercise |
|------|------|----------|
| SC | Standard chow | Sedentary |
| TC | Standard chow | Training |
| SH | High-fat diet (HFD) | Sedentary |
| TH | High-fat diet (HFD) | Training |

Three tissues profiled: **scWAT** (subcutaneous white adipose), **vWAT** (visceral white adipose), **SkM** (skeletal muscle).

Three libraries removed after Step 3 QC: D19-5431, D19-5443, D19-5462.

---

## Bulk RNA-seq Pipeline (independent track)

The bulk pipeline is independent of the scRNA-seq pipeline and validates its findings. A gene that is differentially expressed in bulk and also shows cell-type-specific expression in scRNA-seq is a stronger result than either modality alone.

```
bulk_sample_level_data.txt
      │
      ▼
Filter → protein-coding genes only, deduplicate by MGI symbol
      │
      ▼
Per-tissue DESeq2 loop (scWAT, vWAT, SkM):
      │
      ├── QC assessment
      │     ├── VST normalization (removes mean-variance dependence)
      │     ├── Sample distance heatmap (Ward linkage, Euclidean distance in VST space)
      │     └── PCA on top 500 variable genes (2D colored by pheno_class, 3D by training)
      │
      ├── Differential expression
      │     ├── Design: ~ pheno_class (negative binomial GLM)
      │     ├── Low-count filter: drop genes with total counts < 10
      │     ├── LFC shrinkage: ashr (pulls unreliable low-count estimates toward zero)
      │     └── P-value adjustment: IHW (independent hypothesis weighting by mean expression)
      │
      ├── Three contrasts per tissue
      │     ├── TC vs SC  — exercise effect on healthy diet
      │     ├── SH vs SC  — diet effect (obesity model, sedentary)
      │     └── TH vs SH  — exercise effect in obese animals (key therapeutic contrast)
      │                      requires releveling DESeq2 model to SH as reference
      │
      └── Outputs
            ├── <tissue>_sig_<contrast>_0.05.txt   — significant DE genes (padj < 0.05)
            ├── <tissue>_dist_cluster_protein_coding.pdf
            ├── <tissue>_vsd_pca_pc1-3_protein_coding_top500_training.pdf
            └── dds_<tissue>.pkl
```

---

## Reference Materials

| File | Role |
|------|------|
| `scRNA_analysis_pipeline.ipynb` | Conceptual reference for each step — tool options, method rationale, decision guidance. Not dataset-specific. |
| `README.md` | Pipeline ordering table with one-line descriptions of every file |

---

## Python Package Map

| Stage | R package | Python equivalent |
|-------|-----------|-------------------|
| Ambient RNA removal | SoupX | Manual implementation (Step 3) or CellBender |
| QC + normalization | Seurat | scanpy |
| Doublet detection | DoubletFinder | scrublet |
| tSNE | FIt-SNE | openTSNE |
| Batch correction | Harmony / Seurat rPCA / LIGER | harmonypy / scvi-tools / pyliger |
| Graph clustering | Seurat Leiden | scanpy leiden |
| Density clustering | dbscan | sklearn.cluster.DBSCAN |
| Heatmaps | pheatmap | seaborn.clustermap |
| Bulk DE | DESeq2 + IHW + ashr | pydeseq2 |
| Cell type annotation | SciBet / SingleR | celltypist |
| Cell-cell communication | CellPhoneDB | liana / cellphonedb |
| Gene regulatory networks | SCENIC | pyscenic |
| Copy number variation | inferCNV / CopyKAT | infercnvpy |

---

## File Index

| Step | Notebook | R source | Description |
|------|----------|----------|-------------|
| 1 | `cellranger_metrics.ipynb` | `cellranger_metrics.R` | Collect per-sample CellRanger QC stats |
| 2 | `sample_level_diagnostic_plot.ipynb` | `sample_level_diagnostic_plot.R` | Visualize QC to set filtering thresholds |
| 3 | `sample_level_processing.ipynb` | `sample_level_processing.R` | SoupX → filter → normalize → PCA → DoubletFinder |
| 4 | `sample_level_pseudobulk_clustering.ipynb` | `sample_level_pseudobulk_clustering.R` | Cross-sample pseudobulk sanity check |
| 5 | `sample_integration.ipynb` | `sample_integration.R` + `fancy_tsne.R` | Multi-sample integration → tSNE/UMAP/clustering |
| 6 | `add_metadata.ipynb` | `add_metadata.R` | Annotate cells with phenotype metadata |
| 7 | `analysis_pipeline_scwat_vwat_skm.ipynb` | `analysis_pipeline_scwat_vwat_skm.Rmd` | Downstream analysis (cell types, DEGs, figures) |
| 8 | `8_pseudobulk_differential_expression.ipynb` | — | Pseudobulk DESeq2: obesity + exercise contrasts in vWAT |
| 9 | `9_rescue_and_cell_state_deep_dives.ipynb` | — | Rescue contrast (TH vs SH) + Areg, CP deep dives |
| 10 | `10_diet_ceiling_and_pathway_enrichment.ipynb` | — | Diet ceiling (TH vs TC) + four-contrast classification + GO enrichment |
| — | `scRNA_analysis_pipeline.ipynb` | `scRNA_analysis_pipeline.Rmd` | Pipeline reference documentation |
| — | `bulk_script.ipynb` | `bulk_script.R` | Parallel bulk DESeq2 pipeline (independent) |

---

## NB8 Analysis Plan (`8_pseudobulk_differential_expression.ipynb`)

NB8 implements pseudobulk DESeq2 as the statistically correct alternative to NB7c's Wilcoxon. N=3 mice per condition in vWAT limits power, so fewer DEGs are expected — but what survives is trustworthy.

### Completed

| Section | Status | Key finding |
|---------|--------|-------------|
| Pseudobulk creation | ✓ done | Sum raw counts per mouse × cell state; 6 pseudobulks for SH vs SC |
| EDA (PCA on pseudobulks) | ✓ done | Conditions should separate on PC1; library size check |
| Obesity contrast (SH vs SC) | ✓ done | WAT_IPC: 125 DEGs. Thbs1 confirmed (log2FC=1.7, padj=2e-8). Col3a1 also confirmed. Other ECM genes direction correct but underpowered. |
| Exercise contrast (TC vs SC) | ✓ done | WAT_IPC: 13 DEGs. Nr1d2 confirmed (log2FC=1.3, padj=0.008). Other circadian genes point right direction but N=4 vs 3 insufficient. |
| DEG heatmap (Fig 4G analog) | ✓ done | IPC/CP dominate as expected |
| Wilcoxon vs DESeq2 overlap | ✓ written | Needs NB7c output file to run |

### To add to NB8 (priority order)

**1. Rescue contrast: TH vs SH in WAT_IPC** ← do first

The paper's central therapeutic claim. Already ran this exploratorily — strongest result in the analysis:
- 167 significant DEGs (more than obesity or exercise contrasts alone)
- 5 circadian genes significantly restored: Dbp (log2FC=4.1), Hlf (5.0), Tef (2.3), Nr1d2 (1.9), Per3 (1.7)
- Timp1 significantly reversed: log2FC=-2.7, padj=5e-26
- Arntl (Bmal1) also significantly down (log2FC=-2.9) — a circadian gene going the wrong direction; worth investigating
- Test the paper's "94–95% rescue" claim: of the IPC obesity DEGs (SH vs SC), what fraction flip sign in TH vs SH?

**2. Areg exercise deep dive**

Areg has 30 exercise DEGs — more than WAT_IPC (13) — and is not in the paper's main narrative. The results are surprising:
- Circadian story is present: Dbp (log2FC=4.5, padj=0.03), Tef (1.5), Nr1d2 (1.4), Nr1d1 (1.4) all significant
- **Fndc5 significantly up** (log2FC=1.5, padj=0.008) — Fndc5 is the precursor to irisin, an exercise-induced hormone primarily known as a muscle-secreted myokine. Finding it in adipose stromal cells is unexpected.
- ECM reversal: Col8a1, Postn, Eln, Ccn2 (CTGF — major pro-fibrotic gene) all significantly down with exercise
- Areg = regulatory cells that express amphiregulin; their exercise response has not been described in the paper

**Fndc5/Areg credibility checks — PASSED (2026-05-10)**

Three skeptical hypotheses were tested and all failed to explain the result away:

1. **Muscle contamination?** No — vWAT CP is the #2 Fndc5-expressing cell state across the entire atlas (mean=0.106), higher than Muscle_Fiber (0.071). Fndc5 is expressed broadly, not muscle-specifically. A contamination story predicts near-zero adipose expression; the data shows the opposite.

2. **Single outlier mouse?** No — per-mouse pseudobulk values show perfect rank separation across all 7 animals with no overlap:
   - SC mice: mean log-norm expression 0.020, 0.032, 0.021
   - TC mice: 0.044, 0.063, 0.073, 0.083
   Raw count sums tell the same story: SC mice total 39, 45, 9 Fndc5 counts; TC mice total 49, 78, 64, 69.

3. **Exercise-induced in muscle (contamination direction)?** No — in SkM, Fndc5 goes *down* with exercise in most cell states (Tenocyte_Pappa2+: 0.164→0.082, TSPC: 0.049→0.016). If muscle contamination were driving the Areg TC signal, it would predict a decrease, not an increase.

**Conclusion:** The Fndc5 upregulation in vWAT Areg cells with exercise is a robust signal that survives all three credibility checks. It is not explained by contamination or by a single outlier mouse. Whether this is a genuinely novel finding requires a literature check — Fndc5/irisin in adipose stromal cells may already be documented.

**3. CP obesity deep dive**

CP has 413 obesity DEGs vs IPC's 125 — and it is NOT the same story:
- Top hits: Acta2 (alpha-smooth muscle actin, log2FC=8.4), Tnc (tenascin-C, 5.1), Ltbp2 (latent TGF-beta binding protein, 4.5)
- These are canonical myofibroblast activation markers — IPC's response looks like ECM buildup; CP's response looks like conversion from preadipocyte to fibroblast-like scar-forming cell
- Timp1 and Thbs1 significant in CP obesity too, but Acta2/Tnc/Ltbp2 are not in IPC's top hits at all
- Suggests a differentiation-stage-dependent response: earlier stem cells (IPC) accumulate ECM; more committed cells (CP) trans-differentiate into myofibroblasts

**4. CP exercise deep dive**

CP has 34 exercise DEGs vs WAT_IPC's 13. Circadian genes present but mostly below threshold:
- Tef significant (log2FC=1.5, padj=0.015); Dbp very large effect (5.2) but padj=0.26 — just underpowered
- Nr1d1 significant (1.7, padj=5e-5); Bhlhe41 significant (1.6) — a clock-controlled gene
- Dgat2 significantly up (2.1, padj=0.003) — rate-limiting enzyme for triglyceride synthesis; exercise upregulating fat storage enzyme in committed preadipocytes is counterintuitive and worth a note

### Observations from exploratory analysis

**Tagln vs Timp1 in WAT_IPC obesity (SH vs SC)**

| Gene | log2FC | padj | Per-mouse consistency | Verdict |
|------|--------|------|-----------------------|---------|
| Tagln | 8.2 | 0.018 | Variable — one mouse (D19-5440) drives almost all signal: 468 counts vs 0, 21, 41 in the other two SH mice. SC mice all 0. | **Artifact. Set aside.** |
| Timp1 | 2.7 | 2.7e-21 | Consistent — all 3 SH mice elevated (2006, 1473, 581 counts per mouse); all 3 SC mice lower (360, 216, 55). | **Genuinely interesting.** |

**Timp1** is a direct inhibitor of MMPs (matrix metalloproteinases — enzymes that degrade ECM). When Timp1 is strongly upregulated in IPC cells under HFD, MMP activity is suppressed, ECM that would normally be degraded accumulates, and new ECM from Thbs1/Col genes also accumulates. The fat tissue stiffens and becomes fibrotic. Timp1 is the mechanistic bottleneck in this cascade — not just another ECM gene, but the brake on ECM clearance. The extreme significance (padj=2.7e-21) with clean per-mouse separation makes this the most confident single finding in the obesity contrast. Include prominently in NB8 obesity section.

**Tagln** is a smooth muscle actin-binding protein. The huge log2FC is an N=3 artifact — one mouse drives everything. Do not investigate further.

---

**Mef2c in WAT_IPC exercise (TC vs SC)**

| Gene | log2FC | padj | Per-mouse consistency |
|------|--------|------|-----------------------|
| Mef2c | 4.6 | 0.0003 | Partially consistent — 2 of 4 TC mice drive the signal (D19-5447: 95 counts, D19-5425: 45 counts); D19-5458: 9 counts; D19-5459: 29 counts. All 3 SC mice near-zero (5, 1, 0). |

Mef2c is a transcription factor classically associated with muscle differentiation and cardiac development. It's a downstream target of AMPK and Ca²⁺/calmodulin signaling — both activated by exercise. Finding it as the top exercise-induced gene in fat stem cells raises the hypothesis that exercise is activating muscle-like transcriptional programs in adipose progenitors. The paper foregrounds Nr1d2/circadian genes as the exercise story; Mef2c has a larger effect and better padj but didn't get featured.

**Credibility caveat:** Mef2c's top atlas-wide expressing cell states are B cells (B_Fo mean=1.66 in scWAT), not muscle. This is unexpected for a "muscle TF" — it means Mef2c has broader expression than its classical muscle role suggests. In SkM, it's exercise-responsive in vSMC (+0.45) but not in the FAP or fiber populations where you'd expect it. The per-mouse pattern shows 2 strong TC responders and 2 weak ones, which is less clean than Timp1 or the Fndc5/Areg signal.

**Verdict:** Potentially interesting, needs a literature check before treating as novel. B-cell expression pattern suggests broader function than classical muscle role that may already be documented. Do not lead with this finding; mention as a secondary observation pending literature check.

---

### NB9 findings (2026-05-10)

| Section | Key finding |
|---------|-------------|
| Rescue contrast (TH vs SH) | WAT_IPC: 167 DEGs. Rescue fraction: **96.8%** (121/125 obesity DEGs flip sign — exceeds paper's bulk claim of 94–95%). 5 circadian genes significantly restored. Timp1 significantly reversed (log2FC=-2.7, padj=5e-26). |
| Rescue across cell states | CP: 244 DEGs, Areg: 118, pre_CP: 68 — same cell states dominate as obesity/exercise. Fibroblast: 23 rescue DEGs but only 2 obesity DEGs — exercise-specific effect in obese tissue, not a simple reversal. |
| Areg exercise | ECM reversal is the cleanest signal: Col8a1, Eln, Postn, Thbs1, Ccn2 all significant and negative. Fndc5 confirmed (1.5, padj=0.008). Circadian story: Dbp, Tef, Nr1d1, Nr1d2 all significant. |
| CP obesity | Acta2 (8.4), Tnc (5.1), Ltbp2 (4.5) — myofibroblast activation confirmed. Only 53 genes shared with IPC; 360 CP-unique. |
| CP exercise | Dgat2 (2.1, padj=0.003) confirmed. Also appears in WAT_IPC rescue (3.7, padj=3e-7) — signal is cross-state and cross-contrast, more credible. Dbp large (5.2) but underpowered (padj=0.26). |

**Fibroblast rescue finding:** 23 rescue DEGs (16 down) despite only 2 obesity DEGs. Not a classical rescue — exercise is doing something to fibroblasts in obese tissue that isn't reversing a detectable obesity effect. Likely a sub-threshold obesity activation that falls below N=3 detection but is measurably suppressed by exercise. Top downregulated hits should be checked for fibroblast activation markers (Acta2, Postn, TGF-β targets).

---

## NB10 Plan: Exercise vs Diet — What Can Exercise Not Fix?

**Research question:** TH vs TC isolates the residual effect of HFD in exercising animals. Genes that are still different between TH and TC are things exercise cannot normalize — they require dietary change regardless of fitness. This is the most clinically actionable contrast in the study.

**Contrast:** TH (HFD + exercise) vs TC (low-fat + exercise)
**N:** 3 TH mice vs 4 TC mice

**Why this is distinct from the rescue contrast:** The rescue (TH vs SH) asks "what does exercise do in obese animals?" TH vs TC asks "after exercise, how much residual damage does the diet leave?" A gene that is rescued (TH > SH) but still different from lean exercisers (TH ≠ TC) is one where exercise helps but doesn't fully normalize.

**Analogy:** rescue = treatment response; TH vs TC = treatment ceiling.

### Planned sections

**1. TH vs TC in WAT_IPC** — primary cell state
Run DESeq2 (TH as numerator, TC as reference). Genes significant and positive = HFD still elevating them despite exercise. Genes significant and negative = HFD still suppressing them despite exercise.

**2. Cross-contrast summary for WAT_IPC** — the four-contrast picture
Build a table: for each WAT_IPC gene, show log2FC in all four contrasts (SH vs SC, TC vs SC, TH vs SH, TH vs TC). Classify each obesity DEG as:
- Fully rescued: obese up → rescued down → normalized vs TC
- Partially rescued: obese up → rescued down → still elevated vs TC
- Exercise-resistant: obese up → not rescued → still elevated vs TC

**3. TH vs TC across all vWAT cell states** — heatmap
Same structure as NB8's heatmap. Add TH vs TC column alongside obesity, exercise, and rescue columns for a complete four-contrast picture.

**4. Fibroblast follow-up** — top 16 downregulated rescue DEGs
Check whether fibroblast rescue genes are fibroblast activation markers (Acta2, Postn, TGF-β targets). Run TH vs TC in fibroblasts to see if the sub-threshold obesity effect is visible from the other direction.

### Open questions this will answer

- Are the circadian genes fully normalized by exercise (TH ≈ TC) or just partially rescued (TH between SH and TC)?
- Is Timp1 fully normalized or does some ECM brake persist in TH animals?
- What is the Dgat2 story across all four contrasts — is it an obesity effect, an exercise effect, or both?

---

### Literature check results (2026-05-10, via Google Gemini)

**Literature check completed 2026-05-10 (two passes — corrected after Gemini's circular citation error)**

Gemini's first pass erroneously cited Yang et al. 2022 as prior art. We verified by reading
the paper's CCC section (Figure 5): Yang et al. covered only RANK-RANKL-OPG, MIF-CD74, and
AREG-EGFR. CXCL12 and Fndc5/Areg are absent. Second pass with corrected prompt produced
clean answers on all three questions.

---

**Finding 1: Fndc5/irisin in Areg cells — CONFIRMED NOVEL**
- No paper has reported exercise-induced Fndc5 upregulation in Areg (CD142+) cells in visceral fat.
- Wrann et al. 2013 and follow-up work established Fndc5 in muscle and brain, not CD142+ adipose stroma.
- Broader adipose Fndc5 known (Moreno-Navarrete 2013) but at tissue level in obesity/insulin resistance context, not cell-state specific or exercise-induced.
- Yang et al. 2022 did not report it despite having the data.
- **Verdict: Novel.**

**Finding 2: CXCL12/CCL2 switch — CONFIRMED NOVEL**
- No paper has reported CXCL12 upregulation in visceral Areg/CD142+ cells or committed preadipocytes with exercise.
- Bone marrow precedent exists (LepR+ MSCs, exercise-induced CXCL12) but not in adipose stroma.
- CCL2 decrease with exercise widely known at tissue level, but synchronized CCL2↓/CXCL12↑ in specific adipose stromal cell states has not been described.
- Yang et al. 2022 did not report it despite having the data.
- **Verdict: Novel, though weaker than Finding 1 (CXCL12 failed credibility checks for cell-state specificity and DESeq2 significance). Report as directional/exploratory.**

**Critical new context: Mu et al. 2026**
A 2026 paper (Mu et al., "Irisin ameliorates obesity and insulin resistance via adipose tissue IL-33
and regulatory T cells") established: exogenous irisin → visceral MSCs → IL-33 upregulation →
ST2+ Treg expansion → suppression of metabolic inflammation.

This directly intersects with our Fndc5 finding. We checked Il33 expression in the atlas:

| Cell state | Il33 SC | Il33 TC | Direction |
|------------|---------|---------|-----------|
| Areg | 0.083 | **0.150** | UP with exercise |
| WAT_IPC | 0.069 | 0.076 | flat |
| CP | 0.005 | 0.003 | flat |

**Areg cells increase BOTH Fndc5 AND Il33 with exercise.** Mu et al. showed exogenous irisin
drives Il33 in MSCs. We may be seeing the endogenous version: local Areg-derived irisin
acting autocrinally to drive Il33, which then expands ST2+ Tregs. This extends the NB11
autocrine loop hypothesis with a downstream immune consequence.

Il1rl1 (ST2, the IL-33 receptor) is highest in nILC2 (1.78) and Treg (1.61) — the correct
receiver cells for the Treg expansion story — though expression doesn't change dramatically
with exercise in this dataset.

**This is the most important new finding from the literature check.** The Fndc5→local irisin→
Il33→Treg axis connects our computational finding to a wet-lab validated mechanism from 2026.

**Finding 3: Mef2c in WAT_IPC — Set aside**
B-cell contamination risk, lower priority than Fndc5/CXCL12. Do not pursue.

---

### Open questions

- ~~Does the rescue fraction (TH vs SH) actually hit the paper's claimed 94–95%?~~ **Resolved:** 96.8% in WAT_IPC — exceeds the bulk claim.
- ~~Is Fndc5 expression in Areg cells real biology or a contamination artifact?~~ **Resolved:** signal is real — passes all three credibility checks.
- ~~Is the Fndc5/Areg finding novel?~~ **Resolved:** highly likely novel. Adipose-specific Areg+ source is new; canonical model is muscle-only. Irisin receptor (Itgb5) confirmed expressed in Areg and WAT_IPC, consistent with autocrine loop.
- ~~Does the signal replicate in scWAT or SM FAPs?~~ **Resolved (2026-05-10):** SM FAPs flat (FAP_Areg SC=0.017, TC=0.018, no rank separation) — signal is vWAT-specific, not a general MSC response. scWAT Areg too sparse (4–41 cells/mouse) to confirm or refute — inconclusive, not contradictory. MoTrPAC external validation deprioritized — would only corroborate muscle side, which we already have from SM data.
- ~~Is Mef2c upregulation in WAT_IPC exercise novel?~~ **Partially resolved:** contextually novel (white adipose progenitors, not brown fat), but B-cell contamination risk means it needs validation before claiming. Lower priority than Fndc5.
- Why is Arntl (Bmal1) significantly *down* in the rescue contrast when other circadian genes go up? (Arntl is the positive arm of the clock; Nr1d2/Per3 are the repressive arm — asymmetry may reflect genuine clock phase resetting)
- ~~What is the full Areg exercise signature beyond Fndc5?~~ **Partially resolved by NB12:** CXCL12 is the dominant Areg signal by magnitude; irisin (Fndc5) is real but smaller. See NB12 findings below.
- What are the top 16 downregulated fibroblast rescue DEGs? Are they activation markers? (Needed to distinguish real sub-threshold biology from N=3 artifact)

---

## NB11 Findings: Fndc5/Irisin Autocrine Loop (2026-05-10) ✓ COMPLETE

**Central claim:** vWAT Areg cells are an autonomous exercise-responsive source of irisin, independent of circulating muscle-derived signal.

| Evidence | Key statistic | Interpretation |
|----------|--------------|----------------|
| Fndc5 up in Areg with exercise (vWAT) | log2FC=1.50, padj=0.008 | Only cell state reaching significance |
| Fndc5 down in CP with obesity | log2FC=-1.18, padj=0.003 | Obesity suppresses local circuit |
| Fndc5 rescued in CP (TH vs SH) | log2FC=+1.32, padj=0.004 | Exercise restores circuit in obese animals |
| Itgb5 highly expressed in Areg + WAT_IPC | Mean ~0.52 and ~0.54 | Receptor constitutively present — no barrier to signaling |
| Receptor genes not exercise-regulated | Itgb5, Itgav flat across contrasts | Primed for signaling as soon as ligand appears |
| Fndc5 goes DOWN in SkM with exercise | Tenocyte_Pappa2+: 0.164→0.082 | Cross-tissue asymmetry rules out contamination |
| Per-mouse rank separation in Areg | min(TC)=49 > max(SC)=45 counts | Signal holds in every individual mouse |
| SM FAP states flat with exercise | FAP_Areg SC=0.017, TC=0.018; complete per-mouse overlap | Signal is vWAT Areg-specific — not a general MSC response across tissues |
| scWAT Areg check inconclusive | Group mean direction correct (SC=0.005→TC=0.014) but n=4–41 cells/mouse | Too sparse to confirm or refute; not evidence against |

**Novelty assessment:** Highly likely novel. Gemini deep research found no prior art for Areg-specific exercise-induced Fndc5. Broader adipose Fndc5 is known (Moreno-Navarrete 2013) but not cell-state specific or exercise-induced.

---

## NB12 Findings: Cell-Cell Communication (2026-05-10) ✓ COMPLETE

### Section 1: ADAM proteases
Adam10 and Adam17 constitutively expressed in Areg cells (~0.10–0.14 all conditions). Fndc5 triples TC vs SC (0.024→0.063). Sheddase is a standing capability; exercise flips the substrate. No gap in autocrine loop.

### Section 2: Dominant Areg signal is CXCL12, not irisin
**Biggest surprise from NB12.** Top exercise-amplified interactions from Areg cells are all Cxcl12→Cxcr4, not Fndc5→integrin. Areg increases Cxcl12 signaling to monocytes (+0.35), M1 macrophages (+0.24), T cells, NK cells across the board. This is larger in magnitude than the irisin scores. CXCL12 is a chemokine that directs immune cell positioning — Areg cells are actively repositioning immune populations with exercise.

### Section 3: Irisin paracrine targets are vascular/stromal, not immune
Fndc5→integrin signal increases most in CD4_Memory T cells, vSMC, large vessel, SMC_precursor, pericyte — not macrophages or ILC2s. Receptor driving the paracrine signal is Itgb1 (broadly expressed), not Itgb5 (highly enriched in Areg/WAT_IPC). Somewhat muddies the autocrine specificity story.

### Section 4: ILC2 axis — exercise SUPPRESSES Il33 signaling
Il33→nILC2 signal decreases with exercise (SMC_precursor→nILC2 drops 1.44→0.63 — the largest suppressed interaction in the dataset). Opposite of Gemini's prediction. Biologically coherent: exercise resolves the alarm state rather than amplifying it.

### Section 5: Two clean stromal→immune stories
From stromal→immune pathway heatmap:
- **CXCL12-homing UP**: strongest exercise-amplified pathway across all stromal senders (Areg +1.66, CP +2.34, WAT_IPC +0.96, pre_CP +0.85) — immune cell repositioning
- **Macrophage-recruit (CCL2/CSF1) DOWN**: strongly suppressed — Fibroblast -2.36, WAT_IPC -1.44, Pericyte -0.49 — reduced inflammatory monocyte recruitment

Interpretation: exercise shifts vWAT from CCL2-driven inflammatory macrophage infiltration to CXCL12-driven homeostatic immune positioning.

### Tools ruled out (do not revisit)
- **RNA velocity**: requires spliced/unspliced counts; h5ad only has `layers['counts']`. Would need raw FASTQ re-alignment. Not feasible.
- **CellOracle GRN perturbation**: requires paired ATAC-seq for GRN construction. Dataset is RNA-only. Generic GRN prior would make results entirely assumption-driven.
- **CellPhoneDB**: NB12's manual LR scoring covers the core calculation. With N=3, permutation testing would have near-zero power regardless.

---

## Completed Notebooks

### NB13: Stromal Immune Switch — CXCL12 Up, CCL2 Down ✓ COMPLETE

**Revised framing after credibility checks (2026-05-10):**

Pre-notebook credibility checks ran on CXCL12 analogous to the three checks done for Fndc5.
Results changed the framing significantly:

| Check | Fndc5/Areg result | Cxcl12 result |
|-------|-------------------|---------------|
| Per-mouse rank separation | Clean: min(TC)=49 > max(SC)=45 | Marginal: min(TC)=0.474 ≈ max(SC)=0.468 |
| Cell-state specificity | Areg-specific, only significant cell state | General stromal — CP > Areg, all stromal states move |
| Cross-tissue direction | Opposite in muscle (down) | Same direction in SM (up) — general mesenchymal response |
| DESeq2 significance | padj=0.008, significant | padj=0.27 in Areg, 0.057 in CP — not significant |

**Conclusion from credibility checks:** CXCL12 is not an Areg-specific or vWAT-specific finding.
It is a general mesenchymal/stromal exercise response across tissues. The Fndc5 finding remains
the stronger and more specific of the two.

**Why build NB13 anyway:**
The CXCL12 story is real and directionally consistent even if not cell-state specific. The more
interesting finding is the *combination*: all stromal cells increase CXCL12 (homeostatic immune
positioning) while fibroblasts and WAT_IPC specifically suppress CCL2 (inflammatory recruitment).
That asymmetry — broad CXCL12 up, selective CCL2 down — is a coherent tissue-level immune switch
that is not in the Yang et al. paper and is worth documenting even if it doesn't survive DESeq2.
It also provides context for the Fndc5 finding: Areg cells are doing multiple things with exercise,
and CXCL12 is the larger-magnitude but less specific one.

**Reframed central question:** Does exercise drive a tissue-wide stromal immune repositioning
switch in vWAT — broadly increasing CXCL12 (retention/positioning signal) while selectively
suppressing CCL2 (inflammatory recruitment signal) — and does this switch connect to the
macrophage polarization story from the original paper?

**Questions:**
- All four conditions for Cxcl12 across stromal states — is this lean-exercise specific or also present in rescue (TH vs SH)?
- CCL2 four-condition picture in fibroblasts and WAT_IPC — does obesity increase it, exercise suppress it?
- Cxcr4 receiver landscape — which immune states express the receptor and how does their proportion change?
- Does the CXCL12↑/CCL2↓ switch connect to M1→M2 macrophage polarization (a known exercise effect)?
- Honest framing: what can we claim at mRNA level without DESeq2 significance?

### NB13 Key Findings (2026-05-10)

**1. CCL2 and M1 macrophages track perfectly across all four conditions**
Fibroblast CCL2 and M1 macrophage proportion move in near-perfect lockstep:
- SC: CCL2=0.28, M1=0.25% of cells
- SH (obese): CCL2=0.44, M1=1.75% — obesity turns on fibroblast CCL2, M1s flood in
- TC (lean exercise): CCL2=0.03, M1=0.19% — exercise almost completely suppresses both
- TH (obese exercise): CCL2 partially suppressed, M1 partially reduced

Four-condition dose-response between a stromal ligand and immune cell abundance. Not statistically
tested but pattern across all four conditions is strong. The most visually compelling result in NB13.

**2. Fibroblasts are the dominant CCL2 producers — obesity-activated, not CXCL12 producers**
Fibroblast CXCL12 is near-zero across all conditions. The CCL2/CXCL12 signals come from
different cell states: fibroblasts drive CCL2 (obesity-activated), stromal progenitors (CP, Areg)
drive CXCL12 (exercise-activated). Not the same cells doing both.

**3. Receptor landscape maps cleanly onto the switch**
- CXCR4 highest on Monocyte_Patrolling (1.38 SC, increases in TC) — homeostatic residents
- CCR2 highest on M1 macrophages (1.84) — inflammatory infiltrators
Exercise suppresses CCL2→M1 (inflammatory recruitment) while amplifying CXCL12→patrolling
monocyte (homeostatic retention). Two distinct immune cell populations, two distinct signals.

**Role in overall story:** NB13 provides immune context for the Fndc5 headline. The CCL2/M1
tracking figure is the strongest visual and worth including in any presentation as supporting evidence
that the adipose immune environment is being actively remodeled at the stromal level.

---

### NB14: Fndc5→Il33→Treg Circuit ✓ COMPLETE (2026-05-10)

**Why this jumped to high priority:**
The Mu et al. 2026 paper established: exogenous irisin → visceral MSCs → IL-33 → ST2+ Treg
expansion → suppression of metabolic inflammation. That was shown with injected irisin from
outside the cell.

### NB14 Key Findings (2026-05-10)

**1. Il33 is not a significant DEG in Areg exercise contrast, but trends correctly**
- Areg TC vs SC: log2FC=0.677, padj=0.548 — directionally correct, underpowered
- Fndc5 Areg TC vs SC: log2FC=1.501, padj=0.008 — strong signal (comparison reference)
- Il33 significant in CP: obesity log2FC=2.854, padj=0.032 (up); rescue log2FC=-2.431, padj=7.2e-7 (down)
- Il33 is an active player in the MSC lineage, just not in Areg itself at the pairwise significance level

**2. Fndc5 and Il33 co-upregulate in four-condition pattern in Areg cells**
| Condition | Fndc5 | Il33 |
|-----------|-------|------|
| SC (lean sed.) | 0.024 | 0.083 |
| SH (obese sed.) | 0.014 | 0.122 |
| TC (lean ex.) | 0.063 | 0.150 |
| TH (obese ex.) | 0.021 | 0.116 |
Both genes peak in TC, both suppressed in TH vs TC — same directional logic.

**3. Single-cell co-expression increases 5× with lean exercise**
- SC: 0.3% of Areg cells are double-positive (Fndc5+ and Il33+)
- TC: 1.4% double-positive — 5× increase
- Per-mouse rank separation for Il33 is marginal (one SC mouse overlaps with lowest TC mouse)
  — weaker than Fndc5's clean separation

**4. Treg result is unexpected — honest finding**
| Condition | Treg % of vWAT |
|-----------|----------------|
| SC | 1.58% |
| SH | 3.38% (obesity increases Tregs) |
| TC | 0.78% (lean exercise decreases Tregs) |
| TH | 1.46% |
Tregs increase with obesity, decrease with lean exercise — **opposite of naive Mu et al. prediction**.
Per-mouse pattern is consistent (all TC mice cluster below 1.1%, all SH mice above 2%).
**Interpretation:** Fewer firefighters when there's no fire. TC mice have less adipose inflammation
to suppress, so fewer Tregs needed. The Mu et al. mechanism used obese mice + pharmacological
irisin — likely relevant only in the obesity context. The key Il33 signal in our data is in CP
rescue (padj=7.2e-7), not Areg lean exercise.

**5. ST2 (Il1rl1) is constitutively expressed on Tregs and ILC2s**
- nILC2: mean ~1.7–2.0 across all conditions (slightly elevated in SH)
- Treg: mean ~1.4–2.0 across all conditions
- All stromal states: near zero
- Receptor is available; no barrier to receiving IL-33 signal when produced

**Summary:** The Fndc5→Il33 co-expression in Areg cells is real and consistent, but the Treg
downstream effect is not visible in lean exercise. The strongest Il33 signal is in CP during
obesity and rescue — suggesting the irisin→IL-33 axis may be most relevant when irisin is used
therapeutically in obese animals (matching Mu et al.'s experimental setup), not in lean training.

### NB15: Narrative Document ✓ COMPLETE (2026-05-11)
Standalone narrative notebook for an ML engineer audience. Covers all findings NB11–14,
biological significance, novelty framing, embedded figures, evidence scorecard, honest limits.
Updated after Gemini review to add Section 6c (triple co-expression) and Section 9b (Treg
activation state).

### NB16: Agentic Pipeline Design ✓ COMPLETE (2026-05-10)
Design reference notebook documenting the 6-agent pipeline architecture for automating
this class of re-analysis. Includes working credibility check functions, novelty query
builder, and orchestrator pseudocode. Not executable — design reference only.

---

## Planned Analyses

### Option A: MSC Lineage Specificity Heatmap (LOW EFFORT — 30 min)
**Framing:** Areg is the only vWAT cell state that significantly upregulates Fndc5 with
exercise. But the current evidence is spread across figures. A single dedicated visualization
showing Fndc5 log2FC (and padj markers) across the full Areg → WAT_IPC → pre_CP → CP
differentiation axis makes the "Areg is the specialized irisin-producing stage" claim
explicit rather than implicit.

**What to build:** A focused bar or dot plot of Fndc5 exercise log2FC across the four MSC
lineage stages, with padj stars. Could also include Il33 side-by-side. One panel, ~30 min.

**Narrative gain:** "Only Areg — the most primitive, anti-adipogenic stage — responds with
Fndc5 upregulation. The signal is not present in downstream progenitors (WAT_IPC, pre_CP, CP)
under lean exercise, suggesting Areg cells are specialized irisin producers rather than the
whole MSC lineage coordinating production."

**Where it belongs:** NB14 (new section) or as an additional panel in NB11 Section 3a.

---

### Option B: TF Enrichment via Enrichr/gseapy ✓ COMPLETE (2026-05-11)

**Notebook:** `17_tf_enrichment.ipynb`

**Key findings:**

| Query | Top TF | padj | Interpretation |
|-------|--------|------|----------------|
| Areg up (exercise) | CLOCK | 0.25 (trend) | Fndc5/circadian driven by clock reactivation, not PGC-1α |
| Areg down (exercise) | ESR1 | 2e-6 | Estrogen receptor drives ECM genes exercise suppresses |
| Areg down (exercise) | SMAD4 | 0.009 | TGF-β/SMAD co-drives ECM (Thbs1, Fbn1) |
| CP up (rescue) | SUZ12/PRC2 | 1e-6 | Epigenetic derepression; obesity silences via Polycomb |
| CP down (rescue) | NF-κB1 | 7.6e-8 | Obesity-activated inflammation; exercise suppresses it |
| CP down (rescue) | STAT6 | 5.5e-6 | **Directly includes Il33**; type-2 inflammation driver |

**PGC-1α hypothesis: not confirmed.** The canonical muscle Fndc5 driver does not appear for the Areg exercise program. Instead, Fndc5 is carried along as part of a CLOCK/BMAL1-driven circadian reactivation — the same clock program that drives Nr1d1/Nr1d2/Tef/Dbp. Fndc5 has E-box elements in its promoter, consistent with direct CLOCK regulation.

**Il33 mechanism identified:** In obese CP cells, NF-κB and STAT6 (activated by HFD-driven inflammation) transcriptionally drive Il33. Exercise in obese animals suppresses NF-κB/STAT6 → Il33 falls (padj=7.2e-7). This is the molecular switch connecting exercise to the Mu et al. 2026 IL-33/Treg circuit. The mechanism is primarily relevant in the obese/rescue context, not lean exercise.

**Caveat:** All enrichments from ChIP-seq databases (cell lines, other tissues). Requires ATAC-seq or ChIP-seq in sorted Areg/CP cells to confirm directly.

---

### Option B: TF Enrichment via Enrichr/gseapy (ARCHIVED — see above for results)
**Framing:** We know Fndc5 and Il33 go up in Areg cells with exercise. We don't know *what
is driving them*. TF enrichment analysis on the Areg exercise DEG list can surface candidate
upstream regulators — turning "Fndc5 goes up" into "Exercise activates [TF] in Areg cells,
which drives the Fndc5/Il33 program."

**Approach:** Use `gseapy` to query Enrichr databases (ENCODE ChIP-seq, ChEA, TRANSFAC,
JASPAR) against the Areg TC vs SC significant DEG list. Run the same query on the CP rescue
DEG list (stronger input: Il33 padj=7.2e-7). No new packages beyond gseapy (`pip install
gseapy`).

**Key hypothesis:** PGC-1α (PPARGC1A) is the canonical driver of Fndc5 in muscle. If it
appears as an enriched TF in the Areg exercise signature, that directly connects the fat
and muscle circuits — same regulatory logic, different tissue. Other candidates: NRF2
(oxidative stress response), MEF2 family (exercise-induced in adipose progenitors, NB8).

**CP angle (Option B2):** The CP rescue DEGs include Il33 (padj=7.2e-7) and Fndc5. What
TF drives Il33 suppression in obesity and restoration in exercise in CP cells? If identifiable,
that TF is the mechanistic link between exercise and the Mu et al. IL-33/Treg circuit — a
stronger claim than the co-expression observation alone.

**Output:** A ranked TF list with overlap statistics and a summary of the top 3-5 hits with
biological interpretation. Fits as a new section in NB14 or a standalone NB17.

---

### Option C: Full Regulon Analysis via pySCENIC (HIGH EFFORT — 1 week)
**Framing:** Enrichr (Option B) queries known TF→target relationships from databases.
pySCENIC builds *de novo* regulons from single-cell co-expression patterns, then scores
regulon activity per cell. This is more powerful but requires:
- Download cisTarget databases for mouse mm10 (~15 GB)
- ~2-4 hours compute on 78k vWAT cells
- Interpretation is more complex (regulon activity scores, not simple gene lists)

**When to pursue:** Only if Option B surfaces a candidate TF that's worth validating with
regulon-level evidence. Do Option B first; use pySCENIC to confirm the most interesting hit.

**Fits in:** Standalone NB17 (SCENIC analysis), builds on Option B.

---

### Optional: MoTrPAC Cross-Species Validation
MoTrPAC (NIH exercise atlas, rat endurance training, motrpac-data.org) has bulk RNA-seq
from VENACV (visceral fat) across 1, 2, 4, 8 weeks of training. Could check whether rat
Fndc5 bulk signal goes up in visceral fat with training — cross-species corroboration at
tissue level.

**Key limitation:** Bulk RNA-seq can't confirm the Areg-specific signal (Areg cells are
<1% of total tissue signal). If bulk Fndc5 goes up, it corroborates direction but not
cell-type specificity. If it doesn't show up, that's inconclusive (dilution), not negative.
Honest framing: "Tissue-level corroboration consistent with the mouse Areg signal."

**Estimated effort:** ~2-3 hours. Lower priority than Options A-B since it can't confirm
the core claim.

---

## Presentation Strategy (updated 2026-05-10)

**The headline finding:** Areg cells in vWAT appear to run a self-contained exercise-response
circuit: Fndc5 up → local irisin → Il33 up → Treg expansion. Previously, irisin was thought
to travel from muscle to fat. We show fat may produce it locally, and the downstream consequence
(IL-33/Treg) connects to a wet-lab validated mechanism (Mu et al. 2026).

**How to frame it:**
1. Lead with the Mu et al. 2026 context: irisin→IL-33→Treg is an established wet-lab finding
2. Show our data: Fndc5 and Il33 both go up in Areg cells with exercise (NB11 + NB14)
3. Note what's novel: Mu et al. used exogenous irisin; we see endogenous Fndc5 and Il33
   co-upregulation in the same cell type — consistent with the circuit running locally
4. State the honest limit: this is mRNA, N=3, requires protein-level validation
5. Secondary finding: CCL2↓/CXCL12↑ stromal immune switch (NB13) — directional, not
   statistically confirmed, not in Yang et al.

**What this is NOT:**
- Not a claim that we proved the circuit operates (protein not measured)
- Not a claim that muscle irisin doesn't matter (it might matter too)
- Not a claim that Mu et al. is confirmed (different experimental setup)
