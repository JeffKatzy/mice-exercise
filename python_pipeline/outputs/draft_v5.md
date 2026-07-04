# Exercise co-induces Fndc5 with the circadian clock gene program specifically in visceral fat stem cells

**Jeffrey Katz**

---

## SUMMARY

Yang et al. (*Cell Metabolism*, 2022) established that exercise reactivates circadian clock genes in visceral fat mesenchymal stem cells, identifying this as a headline effect of training. Re-analyzing their published single-cell atlas of 204,883 cells across 51 mice, we find that Fndc5 — the gene encoding the exercise hormone irisin — is co-induced with that same clock program specifically in the most primitive visceral fat progenitors (CD142+ Areg cells), a finding present in their data but not reported. Notably, this co-induction is absent from subcutaneous fat progenitors in mice and undetectable in human subcutaneous fat progenitors, making it depot-specific in a way the CLOCK program alone is not.  

---

## INTRODUCTION

We show that Fndc5 — the gene encoding the exercise hormone irisin — is transcriptionally induced by exercise specifically in the most primitive progenitor cells of visceral fat, and not in subcutaneous fat progenitors in either mice or humans. The induction co-occurs with a circadian clock gene program that was previously identified as a headline exercise effect in these cells, but Fndc5 itself was not reported. The depot-specificity — present in visceral fat, absent from subcutaneous fat — is what makes the observation unusual and worth investigating.

The cells in question are CD142+ Areg cells, first characterized by Schwalie et al. (2018) as gatekeepers at the apex of the visceral fat differentiation hierarchy. Their job is not to become fat cells but to regulate whether the progenitors below them do. The cells directly beneath them — committed preadipocytes — constitutively express Itgb5, the irisin receptor. If Areg cells produce irisin locally during exercise, they would be releasing it directly onto irisin-responsive progenitors at close range — a local signaling architecture that is qualitatively different from circulating irisin arriving from muscle. What that signal would do is not known. But the question is tractable, and this is where it starts.

---

## RESULTS

**Exercise induces Fndc5 specifically in visceral fat Areg cells**

Yang et al. identified differentially expressed genes using cell-level statistical tests, which treat each individual cell as an independent observation. With tens of thousands of cells per condition, this approach is highly sensitive but prone to pseudoreplication — inflating confidence because cells from the same mouse are not truly independent. Applying that same Wilcoxon approach to their published atlas, hundreds of genes reach significance in Areg cells under exercise training; Yang et al. correctly identified the circadian clock genes among them as a headline effect. Applying that approach ourselves, Fndc5 ranks 309th by adjusted p-value in Areg cells — significant, but not highlighted in their analysis. We instead used a pseudobulk approach, summing counts per mouse before testing, which is now the recommended standard for single-cell differential expression. This reduced the significant gene list to 14 — each result supported at the biological replicate level. Running this analysis across all 10,597 expressed genes in vWAT Areg cells, Fndc5 emerged as the thirteenth most significant — not selected in advance, but coming out of an unbiased screen. Among the 14 upregulated genes, four are the canonical CLOCK targets Yang et al. already reported. Fndc5 is the most biologically notable of the remaining ten: it has an established connection to exercise biology, and its receptor Itgb5 is constitutively expressed on the progenitor cells directly below Areg in the hierarchy — an architecture that suggests a local signaling function.

To understand whether the signal was cell-type-specific, we examined Fndc5 expression across all 22 annotated cell types in the vWAT atlas. Only Areg cells showed consistent upregulation under exercise training. Every one of the four trained mice exceeded every one of the three sedentary controls — no exceptions across all pairwise comparisons (logFC=+1.51, padj=0.004; permutation p=0.029). In subcutaneous fat, Fndc5 was near-zero in Areg cells under all four conditions, with no exercise effect. Four of the other 13 upregulated genes are the canonical CLOCK targets — Nr1d1, Tef, Nr1d2, and Dbp — that Yang et al. identified as the core exercise response in these cells (hypergeometric p=3.4×10⁻⁹), confirming our method is recovering the same biology. Fndc5 co-occurs with that program; whether the two are mechanistically linked or simply both triggered by exercise in the same cell is not established by this data.

Two alternative explanations are worth ruling out. The first is technical: could sedentary mice simply have had lower-quality sequencing, producing artificially lower counts? No — sedentary mice had comparable or superior sequencing depth and mitochondrial read rates across all QC metrics (Figure S1). The second is a biological artifact common in single-cell data: ambient RNA. When cells are dissociated from tissue, some rupture and release RNA into the suspension, which gets absorbed by neighboring droplets. If trained mice had more Fndc5-expressing muscle cells nearby, their fat cells could appear to express Fndc5 without actually doing so. But ambient contamination would elevate Fndc5 uniformly across all cell types in trained animals. Instead, elevated Fndc5 was confined entirely to the adipogenic lineage — highest in Areg cells, lower in committed preadipocytes below them, and undetectable in all immune populations such as T cells, B cells, and macrophages (Figure S2). Ambient contamination has no reason to respect cell-type boundaries; a lineage-restricted pattern does. An independent dataset of FACS-sorted CD142+ progenitors, isolated with no surrounding muscle tissue present, confirmed that Fndc5 is genuinely expressed in Areg cells at baseline (GSE128891; Merrick et al., 2019) — showing that the gene is transcriptionally active in these cells even when muscle is entirely absent from the preparation.

In obese exercising mice, all four CLOCK target genes recovered significantly in Areg cells, but Fndc5 did not. This dissociation suggests that Fndc5 induction requires not just an active circadian clock but a second metabolic input — consistent with evidence that AMPK activity, which may be impaired in the obese state, is also required (Guo et al., 2024).

**The obesity-clock suppression link is conserved in humans**

To assess whether this mechanism extends to humans, we analyzed single-nucleus RNA-seq data from 70 donors across lean, obese, and post-bariatric weight-loss groups (GSE295708; Miranda et al., 2025). In Areg-equivalent progenitors from subcutaneous fat, PPARGC1B was significantly suppressed in obesity and restored after weight loss (padj=0.016 and 2.3×10⁻⁷); NR1D2 showed the same pattern after bariatric weight loss (padj=2.9×10⁻¹⁰). FNDC5 was near-zero and not differentially regulated across any comparison in human subcutaneous fat progenitors — mirroring the mouse scWAT null result and consistent with the idea that any FNDC5 response is restricted to visceral fat. Whether exercise induces FNDC5 in human visceral fat progenitors is not addressable from existing public data; no study has combined human visceral fat biopsy with exercise intervention and single-cell profiling. What the human data does establish is that the upstream clock program — suppressed by obesity, recovered by weight loss — is conserved in human fat stem cells, in the same direction as in mice.

---

## DISCUSSION

Yang et al. established circadian clock reactivation in visceral fat MSCs as a headline effect of exercise training. Our contribution is narrower: within that established program, Fndc5 is co-induced specifically in Areg cells, and this co-induction is absent from subcutaneous fat in mice and undetectable in human subcutaneous fat progenitors. The three observations converge — lineage-restricted expression in mice, co-induction with the four canonical clock output genes, and absence from subcutaneous progenitors in both species — on a consistent picture of a fat-autonomous exercise signal concentrated at the top of the visceral fat differentiation hierarchy. Whether the same induction occurs in human visceral fat progenitors is not currently testable from public data.

The depot-specificity is what makes this worth pursuing. Classical thermogenesis via UCP1 operates mainly in subcutaneous fat; visceral fat, the depot most strongly linked to insulin resistance and cardiovascular disease, is largely refractory to it. How exercise affects visceral fat at the cellular level has remained poorly understood. The signal we observe sits at a specific and suggestive location in the tissue architecture: the gatekeeper cells at the very top of the visceral fat hierarchy, whose job is to regulate fat cell formation in the progenitors below them, appear to upregulate Fndc5 during exercise. The committed preadipocytes directly below them constitutively express Itgb5, the irisin receptor. If the mRNA is translated and the protein secreted — which has not been shown — Areg cells would be releasing irisin onto irisin-responsive cells at close range, during exercise, in visceral fat specifically. What that signal would do is genuinely unknown. But the question is tractable: it is one ELISA away from knowing whether the protein is there at all.

In obese exercising mice, all four CLOCK target genes recovered significantly in Areg cells (Dbp padj=1.9×10⁻⁶, Tef padj=1.9×10⁻¹¹, Nr1d2 padj=2.0×10⁻⁹, Nr1d1 padj=0.014), but Fndc5 did not (logFC=+0.65, padj≈1.0). Restoring the circadian program is not sufficient to induce Fndc5 — something present in lean exercise but absent or impaired in obese exercising animals is additionally required. This also clarifies the lean exercise co-induction: Fndc5 and CLOCK genes move together in TC not because one drives the other, but because lean exercise simultaneously satisfies requirements that the obese state separates. Fndc5 is therefore a more selective readout of full exercise adaptation in these cells than CLOCK gene recovery alone.

The immediate next step is direct protein measurement: sort CD142+ Areg cells from trained and sedentary mice and measure FNDC5 protein. Clock dependence can be tested genetically using conditional Bmal1 deletion in F3+ progenitors. The present data is a transcriptional observation pointing toward a mechanism — a fat-autonomous exercise signal, depot-specific and progenitor-restricted, in the tissue where exercise matters most.

---

## METHODS

GSE183288 (Yang et al., 2022); 204,883 cells; log-normalized counts. Areg cells defined by published cell_state_label annotation. Per-mouse expression: log-normalized mean per sample within Areg cells. Permutation test: exact one-sided, all C(7,3)=35 arrangements enumerated; statistic = pairwise inversions. Pseudobulk DESeq2: raw counts summed per mouse (3 SC, 4 TC); genes requiring ≥10 counts in ≥3 samples (10,597 tested); pydeseq2, design `~group`. CLOCK enrichment: hypergeometric test, 16-gene E-box curated target set (Koike et al., 2012). Human validation: GSE295708 (Miranda et al., 2025); Areg-equivalent cells defined by F3>0, PDGFRA>0, PTPRC=0, PECAM1=0; 17 pools; pydeseq2, design `~group + female_only`.

---

## REFERENCES

1. Yang J et al. Single-cell dissection of the obesity-exercise axis in adipose-muscle tissues implies a critical role for mesenchymal stem cells. *Cell Metabolism* 34:1578–1593, 2022.
2. Schwalie PC et al. A stromal cell population that inhibits adipogenesis in mammalian fat depots. *Nature* 559:103–108, 2018.
3. Miranda T et al. A spatially resolved single nucleus atlas of human adipose tissue remodelling in obesity and therapeutic weight loss. *Nature* 627:173–183, 2025.
4. Guo M et al. BMAL1/PGC1α4-FNDC5/irisin axis impacts distinct outcomes of time-of-day resistance exercise. *J Sport Health Sci* 14:100968, 2024.
5. Merrick D et al. Identification of a mesenchymal progenitor cell hierarchy in adipose tissue. *Science* 364:eaav2501, 2019.
6. Koike N et al. Transcriptional architecture and chromatin landscape of the core circadian clock in mammals. *Science* 338:349–354, 2012.
7. Djukic B et al. PKA phosphorylates AMPK to inhibit its activity. *J Biol Chem* 291:17010–17018, 2016.
