"""
Utility for loading DEG ground-truth eval sets across datasets.

Returns a dict keyed by (dataset, tissue, contrast, cell_state) -> list of gene names.

Datasets:
  yang_exercise  — Yang et al. 2022 mouse exercise/obesity scRNA-seq (vWAT, DESeq2 pseudobulk)
  mathys_ad      — Mathys et al. 2019 human AD snRNA-seq (PFC, paper-reported DEGs)
"""

from pathlib import Path
import pandas as pd
from typing import NamedTuple

OUTPUTS_DIR = Path(__file__).parent / "outputs"
MATHYS_DIR  = Path(__file__).parent.parent.parent / "reanalysis_Mathys_2019" / "data"

_YANG_FILES = [
    OUTPUTS_DIR / "degs_vwat_cell_state_deseq2_pseudobulk_sig_only.csv",
    OUTPUTS_DIR / "degs_vwat_cell_state_deseq2_rescue_sig_only.csv",
]

_MATHYS_REPORTED = MATHYS_DIR / "Mathys_reported_DEGs.csv"
_MATHYS_FULL     = MATHYS_DIR / "Mathys_DEGs.csv"
_KANG_SIG        = OUTPUTS_DIR / "degs_kang_pbmc_deseq2_sig_only.csv"
_JERBY_SIG       = OUTPUTS_DIR / "degs_jerby_melanoma_deseq2_sig_only.csv"


class EvalKey(NamedTuple):
    dataset:    str   # "yang_exercise" | "mathys_ad"
    tissue:     str   # "vWAT" | "PFC"
    contrast:   str   # "exercise_TC_vs_SC" | "AD_vs_control"
    cell_state: str   # "Areg" | "Exc" | "Astro" etc.


def _load_yang(files: list[Path] = _YANG_FILES) -> dict[EvalKey, list[str]]:
    frames = [pd.read_csv(f) for f in files]
    df = pd.concat(frames, ignore_index=True)
    result = {}
    for (tissue, contrast, cell_state), group in df.groupby(
        ["tissue", "contrast_label", "cell_state"]
    ):
        key = EvalKey(dataset="yang_exercise", tissue=tissue, contrast=contrast, cell_state=cell_state)
        result[key] = group["gene"].tolist()
    return result


def _load_mathys(path: Path = _MATHYS_REPORTED, padj_threshold: float = 0.05) -> dict[EvalKey, list[str]]:
    df = pd.read_csv(path)

    # Reported DEGs file: use DEGs.Ind.Model flag when present, else filter by adj_p_val
    if "DEGs.Ind.Model" in df.columns:
        sig = df[df["DEGs.Ind.Model"] == True]
    else:
        sig = df[df["adj_p_val"] < padj_threshold]

    result = {}
    for cell_state, group in sig.groupby("celltype"):
        key = EvalKey(dataset="mathys_ad", tissue="PFC", contrast="AD_vs_control", cell_state=cell_state)
        result[key] = group["gene"].tolist()
    return result


def _load_kang(path: Path = _KANG_SIG) -> dict[EvalKey, list[str]]:
    df = pd.read_csv(path)
    result = {}
    for cell_state, group in df.groupby("cell_state"):
        key = EvalKey(dataset="kang_pbmc", tissue="PBMC", contrast="IFNb_stim_vs_ctrl", cell_state=cell_state)
        result[key] = group["gene"].tolist()
    return result


def _load_jerby(path: Path = _JERBY_SIG) -> dict[EvalKey, list[str]]:
    df = pd.read_csv(path)
    result = {}
    for (cell_state, contrast), group in df.groupby(["cell_state", "contrast"]):
        key = EvalKey(dataset="jerby_melanoma", tissue="tumor", contrast=contrast, cell_state=cell_state)
        result[key] = group["gene"].tolist()
    return result


def load_eval_degs(
    include_yang: bool = True,
    include_mathys: bool = True,
    include_kang: bool = True,
    include_jerby: bool = True,
) -> dict[EvalKey, list[str]]:
    """Load all eval sets into a single dict keyed by EvalKey."""
    result = {}
    if include_yang:
        result.update(_load_yang())
    if include_mathys:
        result.update(_load_mathys())
    if include_kang:
        result.update(_load_kang())
    if include_jerby:
        result.update(_load_jerby())
    return result


def summarize(degs: dict[EvalKey, list[str]]) -> pd.DataFrame:
    rows = [
        {"dataset": k.dataset, "tissue": k.tissue, "contrast": k.contrast,
         "cell_state": k.cell_state, "n_degs": len(v)}
        for k, v in degs.items()
    ]
    return (
        pd.DataFrame(rows)
        .sort_values(["dataset", "contrast", "n_degs"], ascending=[True, True, False])
    )


if __name__ == "__main__":
    degs = load_eval_degs()
    print(summarize(degs).to_string(index=False))
    print(f"\nTotal keys: {len(degs)}")
