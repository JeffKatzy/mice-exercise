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
| HFD | High-fat diet — used to induce obesity in mice |
| DEG | Differentially expressed gene — a gene that is significantly more or less active in one condition vs another |
| ECM | Extracellular matrix — the scaffold of proteins between cells; remodeled in obesity |
| UMI | Unique molecular identifier — a short barcode added to each RNA molecule before sequencing to count it exactly once, avoiding PCR duplicates |
| Leiden | A graph-based clustering algorithm; finds "communities" in a cell similarity network |
| DBSCAN | Density-based clustering; finds clusters by local point density, labels outliers as noise |
| tSNE/UMAP | 2D projections for visualizing high-dimensional cell data; preserve local structure |

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
