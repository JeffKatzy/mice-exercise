# Exercise Turns On a Fat Stem Cell Gene — And It's Not Supposed to Be There

**Jeffrey Katz**

*ML engineer's annotated draft — 2026-05-17*

---

## What this is

This is the same paper as `draft_v4.md`, rewritten so I can understand what I actually did. Every section explains the biology problem, the ML/data analogy, what I measured, and what the number means. The scientific claims are identical — only the framing changes.

---

## The Setup: What Kind of Data Is This?

The raw material is a public dataset (GSE183288, Yang et al. 2022) collected by a Harvard/MIT lab. They took 51 mice, split them into four groups — sedentary normal diet (SC), exercised normal diet (TC), sedentary obese (SH), exercised obese (TH) — dissected out fat tissue and skeletal muscle from each mouse, and ran single-cell RNA sequencing on everything.

**What does single-cell RNA sequencing produce?**

Think of it as a feature matrix where:
- Each **row** is a single cell (204,883 rows total)
- Each **column** is a gene (17,341 columns)
- Each **value** is how active that gene is in that cell — essentially how many mRNA transcripts of that gene were captured

This gives you a gene expression profile for every individual cell, not just a tissue average. Instead of "what does fat tissue do under exercise?" (blurry average across all cell types), you can ask "what does *this specific subpopulation* of fat stem cells do?"

The cells are already labeled — the original authors ran clustering (Leiden algorithm, similar to k-means but graph-based) and annotated each cluster as a cell type: macrophages, T cells, fat stem cells, etc. We use their labels rather than re-clustering.

**The four conditions form a 2×2 design:**

```
              Sedentary    Trained
Normal diet:    SC (N=3)    TC (N=4)
High-fat diet:  SH (N=3)    TH (N=4)
```

Our primary contrast is TC vs SC — what does exercise do, independent of obesity?

---

## The Target: What Gene Are We Tracking?

The gene is **Fndc5**. Its protein product is called **irisin** — known as the "exercise hormone" because muscle produces it during exercise, releases it into the bloodstream, and it acts on fat tissue to improve metabolism.

The standard story: **muscle produces irisin → irisin circulates → irisin acts on fat**. Fat is the downstream recipient in this model. Passive.

Our finding: fat stem cells also turn on Fndc5 during exercise. They appear to produce it locally, not just receive it from muscle. Whether this is biologically meaningful depends on whether the mRNA becomes protein and whether that protein acts locally — we haven't shown that yet, but the transcriptional signal is clean.

---

## The Cell Type: What Are Areg Cells?

Fat tissue contains a hierarchy of stem cells, like a differentiation tree:

```
Areg cells (apex — most primitive)
    ↓
WAT_IPC (interstitial progenitor)
    ↓
pre_CP (pre-committed preadipocyte)
    ↓
CP (committed preadipocyte)
    ↓
mature fat cell
```

**Areg cells** (marked by the gene F3/CD142) sit at the very top of this tree. They're unusual: instead of differentiating into fat cells themselves, they *suppress* their neighbors from becoming fat cells. They act as gatekeepers of fat cell formation.

The original Yang et al. paper identified Areg cells as important responders to both obesity and exercise. They found the **CLOCK** transcription factor program gets turned on in fat stem cells broadly by exercise. But they never looked at Fndc5 in Areg cells specifically. That gap is what we exploited.

---

## Finding 1: The Rank Separation

**The question:** Does Fndc5 go up in vWAT Areg cells when mice exercise?

**How we measured it:** For each mouse, we took all the Areg cells from visceral fat (vWAT), and computed the mean log-normalized Fndc5 expression across those cells. One number per mouse.

```
Mouse     Condition    Fndc5
TC_1      Trained      0.044
TC_2      Trained      0.063
TC_3      Trained      0.073
TC_4      Trained      0.083
SC_1      Sedentary    0.020
SC_2      Sedentary    0.032
SC_3      Sedentary    0.021
```

Every trained mouse is above every sedentary mouse. Zero exceptions.

**The statistical test:** With only N=3 vs N=4, we can't run a normal t-test (too few samples, normality assumption meaningless). Instead we use an **exact permutation test**: enumerate every possible way to assign these 7 mice to two groups of size 3 and 4. There are C(7,4) = 35 such arrangements. How many of those 35 arrangements produce a rank separation at least as extreme as what we observed? Answer: 1 out of 35. That's p = 1/35 = 0.029.

This is the minimum achievable p-value at these group sizes. You can't get more consistent than "every trained mouse beats every sedentary mouse" with 7 animals.

**What could make this a false positive?**
- **Bad sequencing quality in SC mice** — lower expression could mean lower read depth, not biology. We checked: SC mice have comparable or *better* sequencing depth than TC mice (SC_1: 2,234 median genes/cell; SC_2: 2,776; SC_3: 1,322). The sedentary mice are not underpowered relative to trained mice.
- **Ambient RNA contamination** — see Finding 2.

---

## Finding 2: Ruling Out Label Leakage (Ambient RNA)

**The problem:** When you dissociate fat tissue into individual cells, some cells rupture and spill their RNA into the suspension. Every droplet picks up a little of this "soup." If the trained mice have higher Fndc5 in muscle cells (which they do), some of that muscle RNA might contaminate fat cells and make it look like fat cells express Fndc5 more.

**This is label leakage** — the label (trained vs sedentary status) is contaminating the features (gene expression) through a technical artifact, not biology.

**How we checked it:** If contamination is the cause, Fndc5 should be elevated uniformly across *all* cell types in trained mice — macrophages, T cells, endothelial cells — not just fat stem cells. Contamination doesn't know which cell type it's landing in.

What we actually see:

```
Cell type         TC mean    SC mean    Interpretation
─────────────────────────────────────────────────────
NK cells          0.000      0.000      ← immune, zero
T cells (CD4/8)   0.000      0.000      ← immune, zero
B cells           0.000      0.000      ← immune, zero
Macrophages       0.000      0.000      ← immune, zero
─────────────────────────────────────────────────────
Areg              0.063      0.024      ← fat stem cell apex
pre_CP            0.019      0.012      ← fat stem, committed
CP                0.154      0.084      ← committed preadipocyte
WAT_IPC           0.003      0.002      ← fat stem, undiff
```

Elevation is restricted entirely to the adipogenic (fat-forming) lineage. The immune cells — which would show contamination if it existed — are flat zero in both conditions. This is the opposite of what contamination looks like.

We also checked an independent orthogonal dataset (Merrick et al. 2019, GSE128891): bulk RNA-seq of FACS-sorted CD142+ progenitors from inguinal adipose tissue shows Fndc5 is detectably expressed at baseline (baseMean = 25.7 normalized counts). This supports that Fndc5 is genuinely present in CD142+ cells and not purely a spillover artifact. Importantly, this dataset is not an exercise intervention, so it supports baseline cell-type expression only, not independent replication of the exercise-induced increase.

---

## Finding 3: The Signal Is Visceral Fat Only

**The dataset has two fat depots:** visceral fat (vWAT, around internal organs) and subcutaneous fat (scWAT, under skin). Same mice, same conditions.

We ran the identical analysis on scWAT Areg cells:
- Most scWAT TC mice: Fndc5 = 0.000
- 3 out of 6 pairwise comparisons show inversions (sedentary mouse higher than trained)
- p = 0.40 — indistinguishable from random chance

We also checked human subcutaneous fat (GSE295708, 70 donors, Miranda et al. 2025): Fndc5 is below the detection threshold entirely (baseMean = 2–3 counts, too low for DESeq2 to test).

So the signal is: vWAT yes, scWAT no (in mice), human SAT no. This pattern is consistent across two species and three independent analyses.

**Why does this matter?** Visceral fat is the metabolically dangerous depot — it drives insulin resistance and cardiovascular disease in obesity. Subcutaneous fat is relatively benign. A response that's specific to visceral fat stem cells is specifically relevant to the fat that actually causes disease.

---

## Finding 4: It's Not Just Fndc5 — There's a Whole Module

**The question:** Is Fndc5 going up by itself, or is it part of a coordinated program?

**Approach:** Run a full differential expression analysis on all 17,341 genes in vWAT Areg cells, comparing TC vs SC at the *mouse* level (not cell level — see the pseudoreplication note below).

**⚠️ The pseudoreplication problem:** The naive approach is to treat each cell as an independent data point — with hundreds of Areg cells per mouse, you'd have N=hundreds and get very low p-values on everything. But those cells are not independent: they're from the same 7 mice. This is like averaging all predictions from one model deployment and treating each request as an independent test of the model. The effective sample size is 7 (the mice), not the number of cells.

**The fix (pseudobulk):** For each mouse, sum the raw gene counts across all its Areg cells into a single "pseudobulk" sample. Now you have 7 samples (3 SC, 4 TC) and run DESeq2 — a proper negative-binomial regression for count data — at the mouse level.

Result: 14 genes significantly upregulated. When we ran the naive cell-level test as a comparison, it yielded 619 genes — a 44× inflation, exactly consistent with the known pseudoreplication problem in single-cell data.

**The CLOCK module:** Among those 14 genes, 4 are canonical outputs of the CLOCK/BMAL1 transcription factor complex:

```
Gene    What it does                           logFC   padj
Nr1d1   Core clock output gene (Rev-erbα)      +1.44   0.0007
Tef     Clock output, metabolic regulator      +1.55   0.005
Nr1d2   Core clock output gene (Rev-erbβ)      +1.37   0.008
Dbp     Strongest known clock output gene      +4.51   0.026
Fndc5   Exercise hormone precursor             +1.51   0.004
```

**What CLOCK/BMAL1 is:** The circadian clock is a gene regulatory circuit that runs on a ~24-hour cycle. The core transcription factors CLOCK and BMAL1 bind to DNA sequences called E-boxes and turn on target genes. These target genes are the "outputs" of the clock — they oscillate daily in healthy tissue. Obesity is known to disrupt circadian gene expression in fat tissue. The original Yang et al. paper showed exercise broadly restores clock genes in fat stem cells.

**The enrichment statistic:** 4 of our 14 significantly upregulated genes are from our curated list of 16 CLOCK targets — out of 10,597 total genes tested. That's a 189-fold enrichment over what you'd expect by chance (hypergeometric p = 3.4×10⁻⁹). The number is extreme, but the more compelling point is just the names: Nr1d1, Nr1d2, Tef, Dbp are *the canonical clock output genes* across tissues and species. These aren't random genes that happen to cluster near a circadian annotation. They are the circadian clock's signature.

Fndc5 is the 6th most significantly upregulated gene in the full dataset, sitting right in the middle of this module.

**A sanity check — why doesn't CLOCK mRNA go up?** CLOCK and BMAL1 proteins are regulated post-translationally — the cell controls their activity through protein phosphorylation and nuclear/cytoplasmic shuttling, not by changing mRNA levels. So flat CLOCK/BMAL1 mRNA with elevated target genes is exactly what you'd expect from increased clock activity. It's like seeing a model making more accurate predictions even though its weights are unchanged — the inputs changed, not the model.

---

## Finding 5: Corroboration From Bulk Data

Single-cell data is noisy — each cell only captures a fraction of its transcripts (dropout problem). So we checked whether the CLOCK module signal is visible in bulk RNA-seq of the same tissue, where you sequence all RNA from the whole tissue homogenate.

Dataset: GSE183239, the companion bulk RNA-seq from the same Yang et al. study. 5 SC vs 5 TC vWAT samples.

```
Gene      SC mean  TC mean   logFC    p
Dbp       2.58     4.97     +2.39    0.008 *
Nr1d2     4.29     5.48     +1.19    0.016 *
Per2      2.48     3.74     +1.25    0.016 *
Per3      3.77     5.18     +1.41    0.032 *
Fndc5     0.00     0.00      0.00    1.0   (undetectable)
```

The CLOCK genes survive bulk dilution. Fndc5 disappears — because Areg cells are less than 1% of total fat tissue, so their Fndc5 signal gets averaged into nothing at the tissue level. This is exactly what you'd expect if the signal is cell-type-specific.

We also checked a Bmal1 knockout microarray dataset (GSE35026, Paschos 2012): Fndc5 in WT adipose doesn't oscillate across circadian times (amplitude 0.17 log2 units across CT0/CT6/CT12/CT18), whereas real clock genes like Nr1d1 have amplitude 2.50. This rules out the interpretation that what we're seeing is just Fndc5 oscillating in sync with the clock. The exercise signal is a training-state effect, not a circadian phase effect.

---

## Finding 6: Human Cross-Validation of the Upstream Mechanism

We can't directly test the visceral fat + exercise question in humans because no public human dataset has visceral fat biopsies + exercise conditions + single-cell resolution. But we can test whether obesity suppresses the same CLOCK genes in human fat stem cells.

Dataset: GSE295708 (Miranda et al. 2025), 70 donors, subcutaneous fat, three groups: lean / obese / post-bariatric weight loss surgery.

We isolated Areg-equivalent cells using the same marker criteria as the mouse analysis (F3>0, PDGFRA>0, immune markers = 0), pooled them into pseudobulk samples, and ran DESeq2.

```
Gene        lean vs obese   padj        wl vs obese    padj
PPARGC1B    +1.02 logFC     0.013 *     +1.84 logFC    4.7e-12 *
NR1D2       +0.68 logFC     0.068       +1.31 logFC    2.4e-13 *
FNDC5       below detection threshold in all groups
```

Two things this tells us:
1. The CLOCK/PPARGC1B program is suppressed by obesity and restored by weight loss in *human* fat progenitors — same direction as the mouse finding, different species, different fat depot.
2. FNDC5 is undetectable in human subcutaneous fat progenitors — consistent with the mouse scWAT null result. The vWAT-specificity of the Fndc5 signal appears to hold across species.

---

## What We're Claiming vs. What We're Not

**What we are claiming:**
- Fndc5 mRNA is elevated in vWAT Areg cells under exercise training, with complete rank separation across 7 mice (p = 0.029, minimum achievable)
- This elevation co-occurs with a coherent set of CLOCK/BMAL1 target genes
- The signal is restricted to the adipogenic lineage (not ambient contamination)
- The signal is absent from subcutaneous fat in both mouse and human
- The CLOCK suppression arm is conserved in human fat progenitors

**What we are NOT claiming:**
- We haven't measured the FNDC5 protein. mRNA going up doesn't guarantee protein goes up. (FNDC5 was undetectable in rat WAT by mass spectrometry, but that's whole-tissue mass spec — it wouldn't detect cell-type-specific expression in <1% of cells anyway.)
- We haven't shown CLOCK directly causes the Fndc5 induction. Clock output genes and Fndc5 going up together is consistent with clock → Fndc5, but it could also be a shared upstream driver. To prove causation you'd need to knock out BMAL1 specifically in F3+ cells and check whether Fndc5 induction disappears. We also looked at a published BMAL1 ChIP-seq dataset (GSE181443) and BMAL1 does not directly bind the Fndc5 gene in inguinal fat — so if CLOCK drives Fndc5 in fat progenitors, it's likely indirect, through downstream nuclear receptors.
- We have N=3–4 mice per condition. The finding is hypothesis-generating, not confirmatory.

---

## What the Experiment That Would Confirm This Looks Like

1. Take trained and sedentary mice
2. Sort vWAT Areg cells by FACS (Lin⁻ CD142⁺ — use the CD142 surface protein as a handle)
3. Measure FNDC5 protein by ELISA on the sorted cells
4. Repeat in F3-Cre × Bmal1-flox conditional knockout mice — if BMAL1 in Areg cells is required for the Fndc5 response, the KO mice should show no induction

This is a tractable wet-lab experiment. The sorting step is standard in adipose biology. The ELISA is routine. The Bmal1-flox mouse line exists.

---

## The Bigger Picture: Why This Might Matter

The metabolically dangerous fat depot (visceral fat) is the one where exercise helps most but the mechanism has been least understood. The classical explanation — exercise burns visceral fat through thermogenic heat generation — barely works in visceral fat; that mechanism is mostly subcutaneous.

This finding points toward a different mechanism: exercise reactivates the circadian clock in visceral fat stem cells, and that reactivation turns on local irisin production in those same cells. Irisin then acts locally (paracrine) on neighboring cells in the most clinically relevant fat depot.

If true, it would suggest:
- Visceral fat has its own exercise response that doesn't depend on receiving irisin from muscle
- The timing of exercise relative to your circadian phase might matter for this specific pathway (there's published data showing this in muscle; untested in vWAT Areg cells)
- This pathway is specifically active under chronic training, not just acute exercise

None of those implications are proven here. This is a computational finding in a public dataset. But it's clean, it's consistent across multiple validation analyses, and it points to a specific experiment that a wet-lab team could run.

---

*Primary dataset: GSE183288 (Yang et al., Cell Metabolism 2022). Human validation: GSE295708 (Miranda et al., Nature 2025). Bulk corroboration: GSE183239. Bmal1 KO reference: GSE35026. CD142+ baseline: GSE128891.*
