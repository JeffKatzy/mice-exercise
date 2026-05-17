# I Dug Through a Dataset from a Top Lab and Found Something They Missed

*A data archeology story about exercise, visceral fat, and a hidden molecular clock*

---

There's a dataset sitting on the NIH's public data servers that contains 204,883 individual cells from 51 mice. It was collected by a team at Harvard and MIT, published in *Cell Metabolism* in 2022, and it cost someone a lot of money and years of work to produce.

The authors found important things. But they didn't find everything.

I found something they missed.

---

## The Setup: Why Visceral Fat Is Different

Before I tell you what I found, you need to understand why visceral fat is the villain in this story.

You have two main types of body fat. **Subcutaneous fat** sits just under your skin — the kind you can pinch. It's metabolically pretty benign. **Visceral fat** wraps around your internal organs. It's the hard belly fat you can't pinch. And it is genuinely dangerous in a way subcutaneous fat is not — it drives insulin resistance, chronic inflammation, cardiovascular disease, and metabolic syndrome. When doctors worry about your waist circumference, they're worried about this stuff.

Here's the frustrating thing: exercise clearly helps visceral fat. But nobody fully understood *how*. The classical explanation — "exercise burns belly fat through thermogenesis" — barely applies to visceral fat. The browning response that generates heat in fat tissue (the UCP1/beige fat pathway you may have heard about) mostly happens in subcutaneous fat. Visceral fat barely does it.

So what is exercise actually doing to visceral fat at the molecular level?

That's what the Yang et al. dataset was built to answer.

---

## The Dataset: A Molecular Atlas of Fat Cells Under Exercise

The original study took mice and put them into four groups:
- **SC**: normal diet, sedentary
- **TC**: normal diet, exercise trained (treadmill, 8 weeks)
- **SH**: high-fat diet, sedentary (the obese group)
- **TH**: high-fat diet, exercise trained

Then they dissected out fat tissue and skeletal muscle, dissociated everything into individual cells, and ran single-cell RNA sequencing on all of them. The result: a gene expression profile for each of 204,883 individual cells, with every cell tagged by tissue, condition, and cell type.

This kind of dataset is enormously powerful. Instead of asking "what does fat tissue do in response to exercise?" — which gives you a blurry average across hundreds of cell types — you can ask "what does *this specific cell type* do?" down to populations of a few thousand cells.

The original authors used this to identify mesenchymal stem cells (MSCs) as key responders to both obesity and exercise. They found that exercise reactivates circadian clock genes in MSCs broadly. Important work.

But they never looked at one specific subpopulation in detail. And that's where things get interesting.

---

## The Hidden Population: Areg Cells

Inside visceral fat tissue, there's a hierarchy of stem cells. At the very top — the most primitive, the ones that sit upstream of everything else — are a population called **Areg cells**. They're marked by a surface protein called CD142 (gene name: F3).

Areg cells are unusual. While most fat stem cells can differentiate into mature fat cells, Areg cells actively *suppress* fat cell formation in their neighbors. They sit at the apex of the differentiation hierarchy and act as gatekeepers. They were first described in a 2018 *Nature* paper from a Swiss group, and the Yang et al. dataset confirmed they exist and respond to both obesity and exercise.

But nobody looked at what specific genes they turn on during exercise training. When I drilled into this population — filtering the dataset down to just the vWAT Areg cells, separating by mouse, separating by condition — a clean signal emerged.

---

## The Finding: Fndc5 Turns On in Visceral Fat Stem Cells During Exercise

**Fndc5** is the gene that encodes irisin — the exercise hormone.

You may have heard of irisin. It's been called the "exercise hormone" because it's produced by muscle during exercise, gets secreted into the bloodstream, and promotes fat browning and metabolic health. There's been some controversy about whether it really circulates in humans (there has, but the current consensus is that it does). Irisin has neuroprotective effects, reduces inflammation in fat tissue, and drives a regulatory T cell response that improves insulin sensitivity.

The standard story: **muscle makes irisin → irisin travels through the blood → irisin acts on fat tissue**.

Fat tissue is the *recipient* in this model. Passive. Downstream.

What I found suggests that's not the whole story — at least for visceral fat stem cells.

When I computed per-mouse mean Fndc5 expression in vWAT Areg cells, here's what the numbers look like:

| Mouse | Condition | Fndc5 expression |
|-------|-----------|-----------------|
| TC_1 | Trained | 0.044 |
| TC_2 | Trained | 0.063 |
| TC_3 | Trained | 0.073 |
| TC_4 | Trained | 0.083 |
| SC_1 | Sedentary | 0.020 |
| SC_2 | Sedentary | 0.032 |
| SC_3 | Sedentary | 0.021 |

Every single trained mouse is above every single sedentary mouse. Zero exceptions. Across all 12 pairwise comparisons, there are zero inversions.

For a dataset with only 3-4 mice per group, that's the cleanest possible signal. The exact statistical test for this (a permutation test that enumerates every possible arrangement of mice to groups) gives p=0.029 — which is the *minimum achievable* p-value at these group sizes. You literally cannot get more consistent than this with N=3 vs N=4.

---

## The First Plot Twist: This Only Happens in Visceral Fat

Here's where it gets more interesting.

The dataset has subcutaneous fat (scWAT) too. Same mice, same conditions. When I ran the same analysis on scWAT Areg cells:

Fndc5 was essentially zero across all conditions. The trained vs sedentary comparison yielded 3 inversions out of 6 pairs — statistically indistinguishable from random chance (p=0.40).

Then I looked at human data. A 2025 *Nature* paper (Miranda et al.) published single-nucleus RNA-seq from subcutaneous fat biopsies of 70 donors — lean, obese, and post-bariatric weight loss surgery. When I isolated the Areg-equivalent cells from that dataset and ran differential expression, FNDC5 was below the detection threshold entirely. Not low — undetectable.

So we have:
- **Mouse visceral fat Areg cells**: clean Fndc5 exercise signal ✓
- **Mouse subcutaneous fat Areg cells**: near-zero, no signal ✗
- **Human subcutaneous fat Areg cells**: below detection ✗

The vWAT-specificity holds across two species and three independent analyses. This isn't an artifact. It's biology.

---

## The Second Plot Twist: It's the Circadian Clock

Why would Fndc5 turn on specifically in visceral fat stem cells during exercise training?

To answer this, I ran a full differential expression analysis comparing trained vs sedentary mice in the vWAT Areg population — testing all 17,341 genes in the dataset. Then I asked: among the genes that go up with exercise, are any of them known targets of the circadian clock transcription factors CLOCK and BMAL1?

I compiled a curated list of 16 genes that are directly controlled by CLOCK/BMAL1 through well-characterized binding sites: Dbp, Nr1d1, Nr1d2, Tef, Per1-3, Cry1-2, and others.

Of the 619 genes significantly upregulated by exercise in vWAT Areg cells, **6 of those 16 CLOCK targets** showed up — a 10.5-fold enrichment over what you'd expect by chance (hypergeometric p=1.2×10⁻⁵).

The single most striking number: **Dbp** — a canonical clock output gene — ranked as the *second most significantly upregulated gene in the entire dataset*, out of 17,341 genes tested. Its fold-change was 4.28 (p_adj = 1.3×10⁻⁸¹).

The picture that emerges: **obesity disrupts the circadian clock in visceral fat stem cells. Exercise restores it. And when the clock turns back on, Fndc5 turns on with it.**

There's supporting evidence for this mechanism from muscle biology. A 2024 paper (Guo et al., *J Sport Health Sci*) demonstrated that BMAL1 — the core clock protein — cooperates with PGC-1α4 to drive Fndc5/irisin expression in skeletal muscle. The mouse Fndc5 promoter has two CLOCK/BMAL1 binding sites. The same machinery that drives irisin in muscle may be operating in visceral fat stem cells, but only when the circadian clock is active.

---

## The Biohacking Angle: What This Actually Means for Your Training

Here's why this matters beyond academic interest.

**1. Visceral fat has its own irisin response — it doesn't just receive irisin from muscle.**

The standard model treats visceral fat as a passive recipient of exercise-induced irisin from muscle. This finding suggests visceral fat stem cells produce irisin locally during training. Local production means it acts on nearby cells in a paracrine fashion — potentially driving effects in the most metabolically dangerous depot that circulating muscle-derived irisin can't reach as effectively.

**2. The circadian clock is the switch.**

This is the finding that connects most directly to the Huberman/Attia world. The upstream event isn't calories burned or some generic exercise signal. It's circadian clock reactivation. Obesity suppresses CLOCK/BMAL1 activity in fat stem cells. Exercise restores it. Dbp, Nr1d1, Nr1d2 — canonical clock genes — are among the strongest exercise-induced genes in the entire fat tissue.

This means **consistent, regular training timing likely matters** for this response. A 2014 study showed circulating irisin peaks at around 9pm in humans. A 2024 study showed that exercising during the late circadian phase — when endogenous BMAL1 activity is highest — significantly enhanced irisin output and downstream metabolic effects compared to early-phase exercise, with the effect abolished in Fndc5-knockout mice. If the same timing dependency applies to visceral fat stem cells (which is testable but not yet shown), then *when* you train may matter as much as *how hard* you train for this specific pathway.

**3. PGC-1β, not PGC-1α.**

A bonus finding: the co-activator that drives the exercise response in fat stem cells is **Ppargc1b** (PGC-1β) — not the PGC-1α you always hear about from muscle biology. In muscle, PGC-1α is already present at baseline and gets activated by exercise through phosphorylation. In visceral fat Areg cells, PGC-1β is essentially absent at baseline and gets transcriptionally induced from near-zero. Different gene, different isoform, different strategy — induced consistently in human subcutaneous fat after resistance exercise too (MoTrPAC dataset, 173 adults, padj=6×10⁻⁵). Fat stem cells are running a distinct molecular exercise program from what muscles run.

---

## What I Actually Did (The Data Work)

For the skeptics: here's what the analysis involved.

I downloaded the original Yang et al. dataset from the NIH GEO database (GSE183288) — the same data the original authors used, converted from their proprietary format to a standard analysis format. No new experiments. No new mice. Just a different set of questions asked of existing data.

For the Fndc5 finding: I computed per-mouse mean log-normalized expression within Areg cells, separated by depot (vWAT vs scWAT), and ran an exact permutation test. I confirmed that QC metrics (sequencing depth, mitochondrial reads) were comparable or higher in sedentary mice — ruling out the possibility that sedentary mice just had worse data quality.

For the ambient RNA check: I verified that the Fndc5 elevation in trained mice is restricted to the adipogenic cell lineage (Areg → committed preadipocytes → progenitors) and is zero in all immune cell populations (NK cells, T cells, B cells, macrophages). Ambient RNA contamination would elevate it everywhere; genuine transcription is cell-type-specific.

For the CLOCK enrichment: I used a curated 16-gene set of direct CLOCK/BMAL1 targets with published binding site evidence, and ran a hypergeometric test against all expressed genes. The enrichment is 10.5-fold with p=1.2×10⁻⁵ — not a fishing expedition.

For the human validation: I analyzed a 2025 *Nature* dataset of 70 human donors (GSE295708), isolated Areg-equivalent cells using the same marker criteria as the mouse analysis, and ran pseudobulk DESeq2. PPARGC1B and NR1D2 are significantly suppressed in obesity and restored by weight loss (padj=0.013 and padj=2.4×10⁻¹³ respectively). FNDC5 is below detection in human subcutaneous fat — which is the expected result given the mouse scWAT data, not a failure to replicate.

---

## What This Doesn't Prove (Being Honest)

This is a computational reanalysis of existing data. The finding is:

- **Hypothesis-generating, not confirmatory.** N=3 mice per condition. p=0.029 is real but the minimum achievable at these group sizes.
- **mRNA only.** We don't know if Fndc5 mRNA upregulation translates to more FNDC5 protein or secreted irisin. That requires sorting Areg cells from trained mice and running an ELISA — a tractable wet-lab experiment that hasn't been done.
- **Correlation, not causation for the CLOCK mechanism.** CLOCK targets co-induce with Fndc5. The Fndc5 promoter has CLOCK binding sites. BMAL1 drives Fndc5 in muscle. But whether CLOCK directly drives Fndc5 in *fat* progenitors requires a ChIP-seq experiment in these cells specifically.

The honest framing: this is a finding that was hiding in a public dataset, is statistically clean, has a plausible mechanism, and makes predictions that a wet-lab team could test. That's what good computational biology is supposed to do.

---

## The Bigger Picture

The original Yang et al. paper cost millions of dollars and years of carefully controlled experiments to produce. The dataset they generated contains more information than any single paper can fully extract. That's normal and expected — big atlases are designed to be mined.

What's interesting about this finding is what it says about visceral fat specifically. We've known for years that exercise is the most effective intervention for visceral fat, and we've struggled to explain *why* at the molecular level. The browning pathway doesn't apply. The explanation for why consistent training specifically shrinks dangerous belly fat has been mechanistically incomplete.

A clock-driven, fat-progenitor-specific irisin response — present in visceral fat stem cells, absent in subcutaneous fat stem cells, triggered by exercise, suppressed by obesity — is a candidate mechanism for something that's been poorly understood for a long time.

The experiment that would prove it is straightforward: sort CD142+ cells from the visceral fat of trained and sedentary mice, measure FNDC5 protein and secreted irisin, knock out BMAL1 conditionally in those cells and see if the response disappears.

Someone should do that experiment.

---

*The analysis code and data are available on GitHub [link]. The original dataset is GSE183288 (Yang et al., Cell Metabolism 2022). Human validation used GSE295708 (Miranda et al., Nature 2025).*
