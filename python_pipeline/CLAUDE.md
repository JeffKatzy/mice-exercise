# Tutorial Pedagogy Guide — mice-exercise Python Pipeline

This file defines the tone, pedagogy, and context rules for all notebooks in
`python_pipeline/`. Use it when updating or adding content.

---

## The reader

An ML engineer learning bio-ML. Strong background in data pipelines, model training,
debugging, and statistics. No assumed biology background. Motivated by understanding
mechanisms, not just running tools.

## Core pedagogy rules

### 1. Problem first, method second
Every major section opens with the problem it solves — not a description of the tool.

BAD:
> SoupX corrects for ambient RNA contamination using three estimation modes...

GOOD:
> When you dissociate tissue to get individual cells, some cells rupture and release
> RNA into the suspension. Every droplet absorbs some of this "soup." Without correction,
> a fat cell will appear to express hemoglobin — which is not real. SoupX estimates how
> much of each cell's observed counts came from the soup and subtracts it.

### 2. Stakes are explicit
For every major step: what specific biological conclusion would be WRONG if you skip it?

Example:
> If you skip doublet removal, two-cell doublets cluster between their parent cell types
> and appear as "rare transitional states." The paper's MSC sub-populations could have
> been entirely artifacts of undetected doublets.

### 3. ML analogies where genuine
Only use analogies that truly map. Label them clearly as analogies.

Good analogies:
- Ambient RNA correction = removing label leakage (contaminating signal from wrong class)
- Doublet detection = finding mislabeled training examples (two samples merged into one)
- Batch correction = domain adaptation (same biology, different technical distributions)
- QC filtering = outlier removal before training
- Pseudobulk = reducing a bag of instances to a bag-level feature vector (MIL)

### 4. Paper-specific numbers are woven in, not appended
Don't say "in this study, 204,883 cells were collected." Say:
> The paper's 204,883 cells came from 39 surviving samples across three tissues —
> three libraries failed and were excluded at this exact step.

### 5. "Why the paper matters" appears near the top of each notebook
One paragraph early on connecting this step to the paper's research question.
The question is: **how does exercise counteract obesity at the cell-type level?**
That question requires resolving individual cell identities, which requires every
cleaning step to work correctly.

---

## Paper context (Yang et al., Cell Metabolism 2022)

**Citation:** Yang J, Vamvini M, Nigro P, et al. Single-cell dissection of the
obesity-exercise axis in adipose-muscle tissues implies a critical role for
mesenchymal stem cells. *Cell Metabolism* 34:1578–1593, 2022.

**Research question:** Exercise training reverses many harmful effects of obesity.
What is happening at the level of individual cell types in fat and muscle tissue?

**Study design:**
- 51 mice, 4 conditions: SC (std chow, sedentary), TC (std chow, trained),
  SH (high-fat diet, sedentary), TH (high-fat diet, trained)
- 3 tissues: scWAT (subcutaneous fat), vWAT (visceral fat), SkM (skeletal muscle)
- 42 samples → 39 after removing 3 failed libraries (D19-5431, D19-5443, D19-5462)
- SVF enrichment used to capture rare stem cells (ASCs in adipose, FAPs in muscle)
  instead of the abundant adipocytes / muscle fibers already well-studied by bulk

**Final atlas:** 204,883 cells, 22 cell types, 42 cell subtypes

**Key finding:** Mesenchymal stem cells (MSCs) are the central responders to both
obesity and exercise — more so than any other cell type. The two main pathways:
- ECM remodeling genes (Thbs1, Sparc) — upregulated by HFD, downregulated by training
- Circadian rhythm genes (Dbp, Tef, Nr1d2, Per3) — upregulated by training in MSCs

**Bulk RNA-seq results:** 1,386 unique DEGs (568 scWAT, 562 vWAT, 256 SkM).
94–95% of adipose DEGs show opposite effects in obesity vs. training ("rescue").

**Three contrasts:**
- Obesity: SH vs. SC (effect of HFD in sedentary animals)
- Training: TC vs. SC (effect of exercise on healthy diet)
- Rescue: TH vs. SH (effect of exercise on already-obese animals) — the key therapeutic contrast

**Human relevance:** Two exercise-regulated genes (DBP and CDKN1A) validated in
independent human cohorts (METSIM study, UK Biobank). MSC exercise signatures
correlate with metabolic traits including BMI and insulin resistance.

---

## Biological vocabulary for non-biologists

| Term | Plain meaning |
|------|--------------|
| scWAT | Subcutaneous white adipose tissue — the fat layer just under your skin |
| vWAT | Visceral white adipose tissue — fat around internal organs; more metabolically harmful |
| SkM | Skeletal muscle |
| SVF | Stromal vascular fraction — everything in fat tissue EXCEPT the mature fat cells; enriched to capture rare stem cells |
| ASC | Adipose stem cell — precursor cells that can become fat cells; key responders to exercise |
| FAP | Fibro-adipogenic progenitor — skeletal muscle stem cell equivalent of ASC |
| MSC | Mesenchymal stem cell — umbrella term for ASCs and FAPs; the paper's key finding |
| IPC | Interstitial progenitor cell — an undifferentiated ASC state |
| CP | Committed preadipocyte — an ASC that has started committing to becoming a fat cell |
| Areg | Adipogenesis-regulatory cell — the most primitive MSC state in fat; marked by F3/CD142; sits at the top of the differentiation hierarchy (Areg → IPC → pre-CP → CP → mature fat cell) |
| HFD | High-fat diet — used to induce obesity in mice |
| DEG | Differentially expressed gene — a gene that is significantly more or less active in one condition vs another |
| ECM | Extracellular matrix — the scaffold of proteins between cells; remodeled in obesity |
| UMI | Unique molecular identifier — a short barcode added to each RNA molecule before sequencing to count it exactly once, avoiding PCR duplicates |
| Leiden | A graph-based clustering algorithm; finds "communities" in a cell similarity network |
| DBSCAN | Density-based clustering; finds clusters by local point density, labels outliers as noise |
| tSNE/UMAP | 2D projections for visualizing high-dimensional cell data; preserve local structure |
| Fndc5/irisin | Fndc5 is the gene; irisin is the protein it encodes after cleavage. Muscle produces irisin during exercise and it circulates in the blood. The Yang et al. finding is that fat stem cells (Areg) also upregulate Fndc5 during exercise — suggesting local irisin production in fat itself, not just import from muscle. |
| AMPK/CaMKII | Exercise sensors in muscle and fat cells. When a cell detects energy stress (AMPK) or calcium flux from contracting muscle (CaMKII), these kinases activate PGC-1α which drives Fndc5 expression. This is the exercise-specific upstream pathway. |
| β3-adrenergic / cAMP | The pathway cold uses to activate fat. Cold → sympathetic nervous system → norepinephrine → β3-AR receptor on fat cells → cAMP → PKA. This is a different upstream route than AMPK/CaMKII. CL316,243 (CL) is a drug that directly activates β3-AR, bypassing the upstream exercise sensors. |
| PGC-1α (Ppargc1a) | A master transcription co-activator for mitochondrial biogenesis and thermogenesis. Sits downstream of both exercise (AMPK/CaMKII) and cold (β3-AR/cAMP). Drives Fndc5, UCP1, and circadian gene expression. |
| Circadian/CLOCK genes | Nr1d1, Nr1d2, Dbp, Tef, Clock, Per1 — genes that run the cell's internal 24-hour clock. HFD disrupts the circadian program in fat cells. Exercise restores it. This is the NB17 finding: CLOCK is the top TF driving the Areg exercise response. |
| UCP1 | Uncoupling protein 1 — the "brown fat" thermogenesis gene. Generates heat by uncoupling the mitochondrial proton gradient. High in BAT, induced in beige fat by cold/exercise. Absent or very low in vWAT progenitors. |
| IL10 | Interleukin-10 — an anti-inflammatory cytokine. In adipose tissue, IL10 from non-thermogenic fat cells suppresses thermogenesis in neighboring cells via IL10R→STAT3 signaling. Feng 2025 showed CLDN5 suppresses IL10 production, enabling thermogenesis. Macrophages also produce IL10 (pro-inflammatory resolution); their IL10 goes the opposite direction from adipocyte IL10 under adrenergic stimulation. |
| Il33/Treg circuit | Il33 is an alarm cytokine produced by stressed CP cells under HFD. It attracts regulatory T cells (Tregs) to fat tissue. Tregs suppress inflammation but also produce IL10, which inhibits thermogenesis. Exercise reduces Il33 in CP cells, unwinding this circuit. This is the NB14 finding. |
| CLDN5/YBX3/IL10 axis | Feng 2025 mechanism: CLDN5 (a tight-junction protein expressed in non-thermogenic fat cells) sequesters YBX3 (a transcription factor) in the cytoplasm. When CLDN5 goes up (cold/exercise), YBX3 stays trapped → IL10 production falls → thermogenic brake is released. When CLDN5 is absent, YBX3 enters nucleus → IL10 rises → thermogenesis suppressed. |

---

## Key causal hypotheses in this project

### Why does exercise upregulate Fndc5 in vWAT Areg cells?

The working hypothesis is a two-step causal chain:

1. **Exercise → CLOCK reactivation in Areg cells.** HFD disrupts circadian rhythms in fat.
   Exercise restores the circadian program (Nr1d1, Nr1d2, Dbp up). NB17 shows CLOCK is the
   top transcription factor driving the Areg exercise response. This is the upstream event.

2. **CLOCK reactivation → Fndc5 upregulation.** CLOCK/BMAL1 directly binds E-box motifs in
   the Fndc5 promoter (established in muscle literature). If the same regulation holds in fat
   progenitors, clock reactivation would drive local Fndc5 expression. This is the mechanistic
   link — not yet proven in fat, which is what makes it interesting.

The β3-agonist experiment (GSE133486) provides a causal test: CL activates cAMP but bypasses
AMPK/CLOCK. Result: circadian genes partially reactivate in APC4 (Areg-like cells), but
Fndc5 stays flat.

**This dissociation is now resolved by Guo et al. 2024** (BMAL1/PGC1α4-FNDC5/irisin axis,
PMC): Fndc5 in muscle requires *both* BMAL1/CLOCK (circadian arm) *and* PGC-1α4 (exercise
energy-sensing arm) co-active simultaneously. β3-agonist activates cAMP/PKA — which wakes
up the circadian clock — but does NOT activate PGC-1α4 (requires actual physical energy
stress/AMPK). So circadian genes reactivate but Fndc5 doesn't follow — exactly what the
two-input model predicts. This makes Fndc5 specifically exercise-gated (not just
adrenergic), which is a strength of the story: only real exercise provides both inputs.

**Gemini literature review (2026-05-15) refined this to a triple-gate model:**
1. **Clock State** — BMAL1:CLOCK active (cAMP/PKA provides this; β3-agonist opens this gate)
2. **Metabolic License** — AMPK active, maintaining Fndc5 promoter poise. Lally et al. 2015:
   AMPK double-knockout dramatically reduces basal Fndc5; AMPK is required but not sufficient alone.
3. **Exercise Key** — PGC-1α4 recruited (physical energy stress / muscle contraction provides
   this; β3-agonist does NOT — it induces PGC-1α1 only)

**Critical finding (Djukic et al. 2016):** PKA actively inhibits AMPK via Ser173
phosphorylation. High-dose β3-agonist may not just fail to activate AMPK — it may actively
suppress it. So β3-agonist provides gate 1, actively blocks gate 2, and provides no gate 3.
Only real exercise provides all three simultaneously.

**MoTrPAC human ASAT corroboration (2026-05-16):** Bulk RNA-seq from 173 sedentary adults
(endurance / resistance / control), subcutaneous fat, acute exercise bouts. Key results:
- **PPARGC1B goes UP** at 3.5-4hr post-resistance exercise: logFC=+0.281, padj=6e-5.
  Endurance trend: +0.171, padj=0.056. Corroborates mouse finding that Ppargc1b is the
  exercise-responsive co-activator in fat, and that resistance exercise drives it more.
- **NR1D1/NR1D2 go DOWN** acutely (NR1D1: −2.2 logFC, padj=1e-48). This is a circadian
  clock reset artifact of acute exercise — genes dip transiently before rebounding. Not
  contradictory to the chronic training upregulation in Yang et al.
- **FNDC5 bulk signal uninformative** — slightly negative but Areg cells are <1% of tissue;
  bulk dilution makes this uninterpretable for cell-type-specific claims.
- **ITGB5 flat** — constitutively expressed, consistent with mouse atlas.

**Ppargc1b finding (2026-05-15):** Ppargc1a4 isoform is unresolvable by scRNA-seq (reads too
short to distinguish isoforms). However, Ppargc1b (PGC-1β) shows the strongest exercise
response of any gene in vWAT Areg cells: SC=0.0010 → TC=0.0065 (6.5× fold change).

**Stromal lineage coordination:** Ppargc1b exercise response spans the entire vWAT MSC
lineage — Areg 6.5×, pre_CP 4.3×, WAT_IPC 3.4×, CP 3.4× — making it a broader exercise
marker across all stromal states. Fndc5 is Areg-specific (only Areg reaches significance);
Ppargc1b is the lineage-wide exercise program marker. CLOCK/BMAL1 drives both as part of a
coordinated program: Fndc5, Ppargc1b, Nr1d1, Nr1d2, Dbp, Tef all move together because
they share E-box binding sites.

**Fat vs muscle mechanistic difference:** Muscle uses Ppargc1a — post-translationally
activated (protein is present at baseline; exercise phosphorylates it into active form). Fat
Areg cells use Ppargc1b, which is near-zero at baseline and transcriptionally induced from
scratch. This is a genuine mechanistic difference: fat stem cells run a different co-activator
program than muscle fibers, which is novel.

**Near-zero per-cell correlation caveat:** Per-cell Spearman r≈0 between Ppargc1b and Fndc5
within Areg cells is due to scRNA-seq dropout (most cells show 0 for both; sparse sampling
makes per-cell correlations uninformative), not mechanistic evidence they're never co-expressed.
The interpretable signal is the group-level four-condition parallel pattern (both peak TC,
both suppressed SH). Correct framing: CLOCK co-induces both in parallel — not that Ppargc1b
directly drives Fndc5 (no evidence for that).

**CLOCK as master regulator:** The exercise program in Areg cells is best described as
CLOCK reactivating a coordinated set of targets: Fndc5, Ppargc1b, Nr1d1, Nr1d2, Dbp, Tef.
Ppargc1b is a secondary novel observation that strengthens the CLOCK-master-regulator
interpretation — both are more compelling together than Fndc5 alone.

### Why would local irisin production in vWAT matter?

vWAT (visceral fat) is metabolically harmful — it drives insulin resistance, inflammation,
systemic disease. Classical thermogenic browning (UCP1-based) barely happens in vWAT.
If Areg cells produce irisin locally, it could act in a paracrine fashion on neighboring
CP cells and preadipocytes — potentially promoting UCP1-independent thermogenesis (P2
futile cycling) or suppressing lipid accumulation. This is why the finding is interesting:
it suggests a local fat-autonomous exercise response in the most clinically relevant depot,
independent of circulating irisin from muscle.

### Human validation findings (2026-05-16)

**GSE214982 — NB19b (omental visceral fat, 9 donors, lean vs obese, Utah):**
- FNDC5 lean=0.0622 vs obese=0.0350 in human omental APC cells (lean 78% higher ✓)
- Direction correct; rank separation imperfect (5 inversions/20 pairs); two female obese donors overlap lean range
- Mouse internal control in same dataset: chow > HFD ✓

**GSE295708 — NB19c v2 (subcutaneous fat, 70 donors, lean/obese/weight_loss, Miranda 2025):**
Pseudobulk DESeq2 on Areg-equivalent cells (F3>0, PDGFRA>0, PTPRC=0, PECAM1=0), 17 pools,
`~ group + female_only` covariate.

| Gene | lean_vs_obese logFC | padj | wl_vs_obese logFC | padj |
|------|---------------------|------|-------------------|------|
| FNDC5 | -0.19 | NaN (below detection) | -0.98 | NaN |
| PPARGC1B | +1.02 | 0.013 * | +1.84 | 4.7e-12 * |
| NR1D1 | +1.76 | 0.166 | +2.30 | 0.001 * |
| NR1D2 | +0.68 | 0.068 | +1.31 | 2.4e-13 * |

**FNDC5 undetectable in human SAT Areg cells** — baseMean=2-3 counts, below DESeq2 threshold.
This is NOT a human replication failure — it is consistent with mouse scWAT biology (see below).

**vWAT-specificity confirmed in mice (checked 2026-05-16):**
Mouse scWAT Areg cells: Fndc5 is near-zero in all conditions (most mice = 0.0000).
Rank separation TC vs SC: 9 inversions/12 pairs — random noise, not a signal.
Mouse vWAT Areg cells: clean rank separation, all 4 TC mice above all 3 SC mice.
Conclusion: Fndc5 exercise response is genuinely vWAT-specific in mice. Human SAT not detecting
it is expected biology, not a technical failure.

**PPARGC1B is vWAT-specific for the exercise response in mice:**
scWAT Areg Ppargc1b: TC mice show 0.000 in 3/4 animals — no exercise response.
vWAT Areg Ppargc1b: 6.2× fold change TC vs SC, clean across all 4 TC animals.
Human SAT PPARGC1B significant difference (lean>obese) reflects obesity-suppression arm,
not exercise induction — consistent with mouse scWAT showing obesity suppresses it there too.

**What the human SAT data does validate:**
CLOCK/circadian program (PPARGC1B, NR1D1/NR1D2) is significantly suppressed by obesity
and restored by weight loss in human SAT Areg-equivalent cells. The upstream mechanism
is conserved; the Fndc5 output is vWAT-specific.

### What makes this novel?

- Irisin/Fndc5 in fat has been described before, but mostly in subcutaneous fat and mostly
  as a response to circulating irisin from muscle.
- Upregulation in vWAT *stem cells* specifically, and at the transcriptional level during
  exercise (not just protein import), has not been reported.
- The CLOCK-driven mechanism in fat progenitors (not muscle) as the upstream driver is new.
- The connection to the Il33/Treg/IL10 thermogenic brake circuit is new.
- **The vWAT-specificity of Fndc5 is confirmed at three levels:** mouse vWAT (clean signal),
  mouse scWAT (near-zero in all conditions), human SAT (undetectable). This is a coherent
  cross-species, cross-depot story — not a replication failure.
- **PPARGC1B and NR1D1/NR1D2 replicate the obesity-suppression arm in human SAT** (DESeq2
  significant), confirming the upstream CLOCK mechanism is conserved even where the Fndc5
  output is depot-restricted.

---

## Structure template for each notebook section

```
### [Section title — name the problem, not the tool]

**The problem:** [1-2 sentences on what goes wrong without this step]

**The approach:** [What we actually do, and why this approach for this data]

**Paper context:** [What Yang et al. did here, what numbers/parameters they used,
what they found as a result of this step]

[Code with inline comments explaining the biological reason for each decision]
```
