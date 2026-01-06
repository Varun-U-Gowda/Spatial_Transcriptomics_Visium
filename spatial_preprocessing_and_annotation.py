#!/usr/bin/env python

"""
Part 1: Cell Type Annotation & Spatial Validation Pipeline (10x Visium)

Author: Varun Umesh Gowda

Purpose
-------
This script performs the first stage of a spatial transcriptomics analysis workflow
for 10x Genomics Visium data. The goal is to generate biologically defensible
cluster- and niche-level annotations by combining:

- Unsupervised clustering (PCA/UMAP + Leiden)
- Canonical marker / mean-expression panels
- Spatial validation (Moran’s I, neighborhood enrichment, spatial gene plots)

Key Outputs
-----------
- Annotated AnnData (.h5ad) with spatial coordinates and clustering results
- Plots: spatial cluster maps, marker gene maps, and UMAP embeddings
- Tables: marker genes per cluster and spatially variable genes

Status
------
Stable / Verified for clustering + marker inspection + TLS-style spatial validation.
Downstream biology (pathways / deconvolution) is handled in Part 2.
"""

from __future__ import annotations

# =============================================================================
# Imports
# =============================================================================

# Standard library
import argparse
import json
import logging
import os
from pathlib import Path

# Third-party
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

# Single-cell / spatial
import scanpy as sc
import squidpy as sq

# ---------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)


# ---------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------
def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Spatial transcriptomics analysis pipeline (Part 1: clustering, DE, SVGs, domains)"
    )
    parser.add_argument(
        "--data_path",
        type=str,
        required=True,
        help="Path to 10x expression matrix (.h5 file)",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default="results/",
        help="Directory to save results",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------
def load_expression_data(h5_path: str):
    """
    Load 10x H5 and auto-use 'spatial/' that sits next to the H5.

    Expected layout:
      /path/to/
        Human_Breast_Cancer_filtered_feature_bc_matrix.h5
        spatial/
          tissue_positions_list.csv
          scalefactors_json.json
          tissue_hires_image.png
          tissue_lowres_image.png
    """
    if not os.path.exists(h5_path):
        raise FileNotFoundError(f"H5 file not found: {h5_path}")

    logging.info(f"Loading 10x expression data from {h5_path}")
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()

    # Derive spatial dir from the H5 path
    base = os.path.dirname(os.path.abspath(h5_path))
    spatial_dir = os.path.join(base, "spatial")
    if not os.path.isdir(spatial_dir):
        raise FileNotFoundError(f"Expected spatial folder at: {spatial_dir}")

    pos_path = os.path.join(spatial_dir, "tissue_positions_list.csv")
    if not os.path.exists(pos_path):
        raise FileNotFoundError(f"Expected positions file at: {pos_path}")

    # Load positions
    positions = pd.read_csv(pos_path, header=None)
    positions = positions.iloc[:, :6]
    positions.columns = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]
    positions["barcode"] = positions["barcode"].astype(str).str.strip()
    positions = positions.set_index("barcode")

    adata.obs_names = adata.obs_names.astype(str).str.strip()

    # First attempt: direct intersection
    common = adata.obs_names.intersection(positions.index)

    # Second attempt: strip trailing "-1" etc.
    if len(common) == 0:
        logging.warning("No overlapping barcodes on first pass; stripping '-#' suffix and retrying.")
        positions.index = positions.index.str.replace(r"-\d+$", "", regex=True)
        adata.obs_names = adata.obs_names.str.replace(r"-\d+$", "", regex=True)
        common = adata.obs_names.intersection(positions.index)

    if len(common) == 0:
        raise ValueError(
            "No overlapping barcodes between H5 and positions. "
            "Confirm that the 'spatial/' folder matches the .h5 file."
        )

    # Subset and attach
    adata = adata[common].copy()
    positions = positions.loc[common]

    adata.obs["in_tissue"] = positions["in_tissue"].astype(int)
    for col in ("array_row", "array_col"):
        if col in positions.columns:
            adata.obs[col] = positions[col].astype(int)

    adata.obsm["spatial"] = positions[
        ["pxl_row_in_fullres", "pxl_col_in_fullres"]
    ].to_numpy(dtype=float)

    # Images and scalefactors
    sf_path = os.path.join(spatial_dir, "scalefactors_json.json")
    scalefactors = {}
    if os.path.exists(sf_path):
        try:
            with open(sf_path, "r") as f:
                scalefactors = json.load(f)
        except Exception as e:
            logging.warning(f"Could not load scalefactors_json.json: {e}")

    images = {}
    for name in ("tissue_hires_image.png", "tissue_lowres_image.png"):
        p = os.path.join(spatial_dir, name)
        if os.path.exists(p):
            try:
                images["hires" if "hires" in name else "lowres"] = np.array(Image.open(p))
            except Exception as e:
                logging.warning(f"Could not load image {p}: {e}")

    adata.uns.setdefault("spatial", {})
    # use a generic sample key
    adata.uns["spatial"]["sample"] = {
        "images": images,
        "scalefactors": scalefactors,
        "metadata": {
            "source_image_path": os.path.join(spatial_dir, "tissue_hires_image.png")
            if os.path.exists(os.path.join(spatial_dir, "tissue_hires_image.png"))
            else None
        },
    }

    logging.info(f"Using spatial directory: {spatial_dir}")
    logging.info(f"Loaded expression + spatial metadata for {adata.n_obs} spots")
    logging.info(f"In-tissue spots: {(adata.obs['in_tissue'] == 1).sum()}")
    logging.info(f"Genes detected: {adata.n_vars}")

    return adata


# ---------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------
def preprocess(adata):
    """
    Filter to in-tissue spots, add mito QC, light QC.
    """
    # Keep only in-tissue spots
    if "in_tissue" in adata.obs:
        adata = adata[adata.obs["in_tissue"] == 1].copy()
    else:
        logging.warning("'in_tissue' not found in adata.obs; proceeding without in-tissue filtering.")

    # Mito genes (human: 'MT-')
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    # Basic QC metrics incl. mito%
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Light QC thresholds (tune if needed)
    if "n_genes_by_counts" in adata.obs:
        adata = adata[adata.obs["n_genes_by_counts"] >= 500].copy()
    else:
        logging.warning("'n_genes_by_counts' not in adata.obs; skipping gene-count filter.")

    if "pct_counts_mt" in adata.obs:
        adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
    else:
        logging.warning("'pct_counts_mt' not in adata.obs; skipping mito% filter.")

    logging.info(f"After QC: {adata.n_obs} spots, {adata.n_vars} genes")
    return adata


def normalization(adata):
    """Normalize and log-transform."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    logging.info("Normalization complete")
    return adata


def select_hvgs(adata):
    """Select highly variable genes."""
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    if "highly_variable" not in adata.var:
        raise RuntimeError("Failed to compute highly_variable genes; check input data.")
    hvgs = int(adata.var["highly_variable"].sum())
    logging.info(f"Selected {hvgs} highly variable genes")
    return adata


def run_pca_umap_clustering(adata):
    """PCA, neighbors, UMAP, and Leiden clustering."""
    sc.tl.pca(adata, n_comps=30, use_highly_variable=True, svd_solver="arpack")
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=15)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    n_clusters = adata.obs["leiden"].nunique()
    logging.info(f"Clustering complete: {n_clusters} clusters identified")
    return adata


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
def plot_spatial_maps(adata, out_dir: str):
    """
    Save cluster map and QC maps on the tissue.
    """
    os.makedirs(out_dir, exist_ok=True)

    # Cluster map
    sc.pl.spatial(adata, color="leiden", spot_size=None, show=False)
    plt.savefig(os.path.join(out_dir, "spatial_clusters.png"), bbox_inches="tight", dpi=300)
    plt.close()

    # QC maps: counts / genes / mito%
    for key in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        if key in adata.obs:
            sc.pl.spatial(adata, color=key, spot_size=None, show=False)
            plt.savefig(os.path.join(out_dir, f"spatial_{key}.png"), bbox_inches="tight", dpi=300)
            plt.close()

    logging.info(f"Spatial maps saved to {out_dir}")


def plot_umap_leiden(adata, out_dir: str):
    """
    Save UMAP with Leiden clusters only (no cell types yet in Part 1).
    """
    os.makedirs(out_dir, exist_ok=True)
    sc.pl.umap(adata, color=["leiden"], show=False)
    plt.savefig(os.path.join(out_dir, "umap_leiden.png"), bbox_inches="tight", dpi=300)
    plt.close()
    logging.info("UMAP (Leiden) saved.")


# ---------------------------------------------------------------------
# Marker genes (DE per cluster)
# ---------------------------------------------------------------------
def find_and_save_markers(adata, out_dir: str):
    """Identify upregulated marker genes per cluster and save filtered table."""
    os.makedirs(out_dir, exist_ok=True)
    logging.info("Finding marker genes (DE per cluster)...")

    if "leiden" not in adata.obs:
        raise RuntimeError("No 'leiden' column found in adata.obs. Run clustering first.")

    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    marker_df = sc.get.rank_genes_groups_df(adata, group=None)

    required_cols = {"group", "names", "logfoldchanges", "pvals_adj"}
    missing_cols = required_cols - set(marker_df.columns)
    if missing_cols:
        raise RuntimeError(f"rank_genes_groups results missing columns: {missing_cols}")

    # Filter: upregulated only
    marker_df_filtered = marker_df[
        (marker_df["pvals_adj"] < 0.05) & (marker_df["logfoldchanges"] > 1)
    ].copy()

    if marker_df_filtered.empty:
        logging.warning("No upregulated markers passed the filters (padj < 0.05, logFC > 1). Saving full table.")
        marker_df.to_csv(
            os.path.join(out_dir, "cluster_markers_all.csv"),
            index=False,
        )
        return marker_df

    # Top 20 per cluster
    all_markers = []
    for cluster_id in sorted(marker_df_filtered["group"].unique()):
        cluster_markers = marker_df_filtered[marker_df_filtered["group"] == cluster_id]
        top20 = cluster_markers.nlargest(20, "logfoldchanges")
        all_markers.append(top20)

    final_markers = pd.concat(all_markers, axis=0)

    # Round numeric columns if present
    for col in ["logfoldchanges", "pvals", "pvals_adj", "scores"]:
        if col in final_markers.columns:
            final_markers[col] = final_markers[col].round(4)

    out_path = os.path.join(out_dir, "cluster_markers_upregulated.csv")
    final_markers.to_csv(out_path, index=False)
    logging.info(f"Marker genes saved to {out_path}")

    return final_markers


# ---------------------------------------------------------------------
# Spatially variable genes (Moran's I)
# ---------------------------------------------------------------------
def find_spatially_variable_genes(adata, out_dir: str):
    """
    Find genes with spatial expression patterns using Moran's I.
    These genes vary across tissue space, not just between clusters.
    """
    os.makedirs(out_dir, exist_ok=True)
    logging.info("Finding spatially variable genes (SVGs) with Moran's I...")

    # Ensure spatial neighbors are computed
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)

    # Calculate Moran's I for spatial autocorrelation
    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        n_perms=100,
        n_jobs=4,
    )

    if "moranI" not in adata.uns:
        logging.warning("Moran's I results not found in adata.uns['moranI']. Skipping SVG export.")
        return adata

    svg_results = adata.uns["moranI"].copy()

    # Choose a p/FDR column robustly
    q_candidates = ["pval_norm_fdr_bh", "pval_z_sim_fdr_bh", "pval_sim_fdr_bh"]
    p_candidates = ["pval_norm", "pval_z_sim", "pval_sim"]

    q_col = next((c for c in q_candidates if c in svg_results.columns), None)
    p_col = next((c for c in p_candidates if c in svg_results.columns), None)

    if q_col:
        svg_results = svg_results.sort_values(q_col, ascending=True)
        svg_results_sig = svg_results[svg_results[q_col] < 0.05].copy()
        svg_results_sig.rename(columns={q_col: "qval"}, inplace=True)
        if p_col:
            svg_results_sig.rename(columns={p_col: "pval"}, inplace=True)
    elif p_col:
        logging.warning("FDR column not found in Moran results; using raw p-values instead.")
        svg_results = svg_results.sort_values(p_col, ascending=True)
        svg_results_sig = svg_results[svg_results[p_col] < 0.05].copy()
        svg_results_sig.rename(columns={p_col: "pval"}, inplace=True)
    else:
        logging.warning(
            "No p-value columns found in Moran results. "
            "Falling back to sorting by |I| and exporting top genes."
        )
        if "I" not in svg_results.columns:
            logging.error("Moran results missing 'I' column. Skipping SVG export.")
            return adata
        svg_results_sig = svg_results.reindex(svg_results["I"].abs().sort_values(ascending=False).index).copy()

    if svg_results_sig.empty:
        logging.warning("No significant spatially variable genes found at FDR/p < 0.05.")
    else:
        out_csv = os.path.join(out_dir, "spatially_variable_genes_moran.csv")
        svg_results_sig.to_csv(out_csv)
        logging.info(f"Significant Moran SVGs saved to {out_csv}")

    # Get top 20 for printing
    if "I" in svg_results_sig.columns:
        top_svgs = svg_results_sig.sort_values("I", ascending=False).head(20)
        logging.info("Top spatially variable genes (by Moran's I):")
        logging.info(top_svgs[["I"] + [c for c in ["qval", "pval"] if c in top_svgs.columns]].to_string())
    else:
        top_svgs = svg_results_sig.head(20)

    # Visualize top 6 SVGs on tissue
    if not top_svgs.empty:
        top_genes = [g for g in top_svgs.index if g in adata.var_names]
        top6_genes = top_genes[:6]

        if top6_genes:
            fig, axes = plt.subplots(2, 3, figsize=(15, 10))
            axes = axes.flatten()

            for idx, gene in enumerate(top6_genes):
                title_extra = ""
                if "I" in svg_results_sig.columns:
                    title_extra = f" (I={svg_results_sig.loc[gene, 'I']:.3f})"
                sc.pl.spatial(
                    adata,
                    color=gene,
                    ax=axes[idx],
                    show=False,
                    title=f"{gene}{title_extra}",
                )

            plt.tight_layout()
            out_png = os.path.join(out_dir, "top_svgs_spatial.png")
            plt.savefig(out_png, bbox_inches="tight", dpi=300)
            plt.close()
            logging.info(f"Top SVG spatial plots saved to {out_png}")

    return adata


# ---------------------------------------------------------------------
# Spatial domains
# ---------------------------------------------------------------------
def identify_spatial_domains(adata, out_dir: str):
    """
    Identify spatial domains/niches - regions with similar expression AND proximity.
    """
    os.makedirs(out_dir, exist_ok=True)
    logging.info("Identifying spatial domains using spatial neighbors...")

    # Ensure spatial neighbors are computed
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)

    # Cluster using spatial constraints (higher weight on spatial neighbors)
    sc.tl.leiden(
        adata,
        resolution=0.3,
        key_added="spatial_domain",
        obsp="spatial_connectivities",
    )

    n_domains = adata.obs["spatial_domain"].nunique()
    logging.info(f"Identified {n_domains} spatial domains")

    # Visualize spatial domains
    sc.pl.spatial(
        adata,
        color="spatial_domain",
        spot_size=None,
        show=False,
        title="Spatial Domains",
    )
    out_png = os.path.join(out_dir, "spatial_domains.png")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close()

    logging.info(f"Spatial domains plot saved to {out_png}")
    return adata


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    """Run Part 1 end-to-end: load Visium data, cluster, inspect markers, and run spatial validation.

    Notes
    -----
    - Visium spots can contain mixed cell populations; interpret markers in a niche-aware manner.
    - This script produces an annotated AnnData for downstream analyses (Part 2).
    """

    args = parse_args()
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    try:
        # 1) Load data + spatial metadata
        adata = load_expression_data(args.data_path)

        # 2) Preprocess & normalize
        adata = preprocess(adata)
        adata = normalization(adata)
        adata = select_hvgs(adata)

        # 3) PCA / neighbors / UMAP / Leiden
        adata = run_pca_umap_clustering(adata)

        # 4) Plots for clusters & QC
        plot_spatial_maps(adata, out_dir)
        plot_umap_leiden(adata, out_dir)

        # 5) DE markers per cluster
        find_and_save_markers(adata, out_dir)

        # 6) Spatially variable genes (Moran SVGs)
        find_spatially_variable_genes(adata, out_dir)

        # 7) Spatial domains based on spatial connectivities
        identify_spatial_domains(adata, out_dir)

        # 8) Save processed AnnData
        out_h5ad = os.path.join(out_dir, "adata_spatial_processed.h5ad")
        adata.write(out_h5ad)
        logging.info(f"Final processed AnnData saved to {out_h5ad}")

    except Exception as e:
        logging.exception("Spatial pipeline failed.")
        # Re-raise so that in a workflow manager this still fails visibly
        raise


if __name__ == "__main__":
    main()
