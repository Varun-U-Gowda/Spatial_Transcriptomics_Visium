#!/usr/bin/env python3
"""
TLS-like validation + (optional) spatially variable genes (Moran's I) for Visium (10x).

What it does
------------
- Loads Visium AnnData (spots + histology coordinates)
- Computes immune/TLS module scores (Scanpy score_genes)
- Creates spatial "tissue vs dots" plots:
    Left  = tissue only (Scanpy)
    Right = dots-only scatter (matplotlib) for guaranteed visibility
- Runs targeted Moran’s I on TLS-related genes
- Runs neighborhood enrichment (cluster adjacency) + saves plot
- Exports summary tables

Notes
-----
- Designed for macOS/Python 3.11 (multiprocessing safe via __main__ guard)
- Older Squidpy versions store Moran results in adata.uns["moranI"]
- Visium spots are mixtures; these plots validate *niche-level* spatial organization
"""

import os
import argparse
import logging
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description=(
            "Spatial plotting/validation utilities for 10x Visium: "
            "module scoring + tissue-vs-dots plots + Moran's I + neighborhood enrichment."
        )
    )
    p.add_argument(
        "--adata_path",
        default="results_auto/adata_spatial_processed.h5ad",
        help='Path to input AnnData (.h5ad) with .obsm["spatial"].',
    )
    p.add_argument(
        "--outdir",
        default="results_auto/tls_validation",
        help="Output directory for plots and CSVs.",
    )
    p.add_argument(
        "--cluster_key",
        default="leiden",
        help="Column in adata.obs containing cluster labels.",
    )
    p.add_argument(
        "--target_cluster",
        default="0",
        help="Cluster label to highlight for cluster-specific summaries/plots.",
    )
    p.add_argument(
        "--cluster_markers_csv",
        default=None,
        help="Optional CSV of cluster-wise upregulated markers to summarize (optional).",
    )

    # Plot controls
    p.add_argument("--dots_marker_size", type=float, default=6.0, help="Marker size for dot plots.")
    p.add_argument("--vmax_gene", default="p99", help="vmax for gene plots (e.g., p99, p98, 3.0).")
    p.add_argument("--vmin_gene", type=float, default=0.0, help="vmin for gene plots.")
    p.add_argument("--vmax_score", type=float, default=0.5, help="vmax for score plots.")
    p.add_argument("--vmin_score", type=float, default=-0.5, help="vmin for score plots.")
    p.add_argument("--dpi", type=int, default=170, help="Figure DPI for saved plots.")

    # Neighborhood enrichment controls
    p.add_argument("--nhood_n_perms", type=int, default=1000, help="Permutations for neighborhood enrichment.")
    p.add_argument("--nhood_n_jobs", type=int, default=1, help="Parallel jobs for neighborhood enrichment.")

    # Optional genome-wide SVG discovery (slow)
    p.add_argument(
        "--run_svg",
        action="store_true",
        help="If set, run Moran's I on HVGs to discover spatially variable genes (slow).",
    )
    p.add_argument("--hvg_n", type=int, default=2000, help="Number of HVGs if --run_svg is set.")
    p.add_argument("--svg_top_n", type=int, default=200, help="Top N SVGs to export if --run_svg is set.")

    return p.parse_args()


def setup_logging():
    """Configure logging."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")


def _setup_scanpy(dpi=170):
    """Global Scanpy plotting defaults (safe for CLI runs)."""
    sc.set_figure_params(dpi=dpi, frameon=False)
    warnings.filterwarnings("ignore", category=UserWarning)


# =============================================================================
# Helpers
# =============================================================================

def present_genes(adata, genes):
    """Return list of genes that are present in adata.var_names (preserves input order)."""
    avail = set(adata.var_names)
    return [g for g in genes if g in avail]


def _get_vals_1d(adata, key):
    """
    Return 1D float array for a gene (adata.var_names) or obs column (adata.obs).

    Returns
    -------
    (vals, src) where src is "gene" or "obs".
    """
    if key in adata.obs.columns:
        vals = adata.obs[key].to_numpy(dtype=float)
        return vals, "obs"
    if key in adata.var_names:
        x = adata[:, key].X
        if hasattr(x, "toarray"):
            x = x.toarray()
        vals = np.ravel(x).astype(float)
        return vals, "gene"
    raise ValueError(f"[ERROR] '{key}' not found in adata.obs or adata.var_names")


def _resolve_vmax(vals, vmax):
    """Support vmax='p99' style percentile strings, otherwise numeric."""
    if isinstance(vmax, str) and vmax.lower().startswith("p"):
        p = float(vmax[1:])
        return float(np.nanpercentile(vals, p))
    return float(vmax)


def count_nonzero(adata, gene):
    """Count nonzero expression entries for a gene across spots (or return None if gene absent)."""
    if gene not in adata.var_names:
        return None
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return int((x > 0).sum())


def spatial_tissue_vs_dots(
    adata,
    color,
    outdir,
    fname_prefix="spatial_gene_localization",
    dots_marker_size=12,
    vmin=0.0,
    vmax="p99",
    dpi=170,
):
    """
    Save a side-by-side spatial plot:

      Left:  tissue only (Scanpy)
      Right: dots-only (matplotlib scatter) for guaranteed visibility + explicit scale.

    Output:
      {outdir}/{fname_prefix}_{color}.png
    """
    coords = adata.obsm["spatial"]
    x = coords[:, 0]
    y = coords[:, 1]

    vals, src = _get_vals_1d(adata, color)
    vals = np.where(np.isfinite(vals), vals, np.nan)

    vmax_val = _resolve_vmax(vals, vmax)
    vmin_val = float(vmin)

    fig, (A, B) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: tissue only
    sc.pl.spatial(adata, color=None, ax=A, show=False)
    A.set_title("Tissue")

    # Right: dots only
    sc_ = B.scatter(x, y, c=vals, s=dots_marker_size, vmin=vmin_val, vmax=vmax_val)
    B.set_title(f"{color} ({src})")
    B.set_aspect("equal")
    B.invert_yaxis()  # Visium coordinate orientation
    B.axis("off")
    plt.colorbar(sc_, ax=B, fraction=0.046, pad=0.04)

    plt.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"{fname_prefix}_{color}.png")
    fig.savefig(outpath, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return outpath


def spatial_tissue_vs_clusters(
    adata,
    cluster_key,
    outdir,
    fname,
    target_cluster=None,
    dots_marker_size=12,
    dpi=170,
):
    """
    Save a side-by-side cluster plot:

      Left:  tissue only
      Right: dots-only colored by cluster (categorical)

    If target_cluster is provided, all other clusters are shown in gray.
    """
    coords = adata.obsm["spatial"]
    x = coords[:, 0]
    y = coords[:, 1]

    clusters = adata.obs[cluster_key].astype(str).values
    unique_clusters = np.unique(clusters)

    cluster_to_int = {c: i for i, c in enumerate(unique_clusters)}
    ints = np.array([cluster_to_int[c] for c in clusters])

    fig, (A, B) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: tissue only
    sc.pl.spatial(adata, color=None, ax=A, show=False)
    A.set_title("Tissue")

    # Right: dots-only
    if target_cluster is None:
        cmap = plt.cm.get_cmap("tab20", len(unique_clusters))
        B.scatter(x, y, c=ints, s=dots_marker_size, cmap=cmap)

        handles = [
            plt.Line2D(
                [0], [0],
                marker="o",
                color="w",
                label=c,
                markerfacecolor=cmap(cluster_to_int[c]),
                markersize=6,
            )
            for c in unique_clusters
        ]
        B.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc="upper left")
        B.set_title(f"{cluster_key} (all clusters)")
    else:
        colors = np.array(["lightgray"] * len(clusters), dtype=object)
        colors[clusters == str(target_cluster)] = "red"
        B.scatter(x, y, c=colors, s=dots_marker_size)
        B.set_title(f"{cluster_key} (cluster {target_cluster})")

    B.set_aspect("equal")
    B.invert_yaxis()
    B.axis("off")

    plt.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, fname)
    fig.savefig(outpath, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return outpath


def export_tls_de_check(cluster_markers_csv, outdir, target_cluster, tls_genes):
    """
    Optional: if you have cluster-wise DE markers (rank_genes_groups on leiden),
    export the subset of TLS genes for the target cluster.
    """
    if not cluster_markers_csv or not os.path.exists(cluster_markers_csv):
        logging.info("cluster_markers_upregulated.csv not provided/found; skipping DE check.")
        return

    logging.info("Reading cluster markers CSV (optional DE sanity-check)...")
    cm = pd.read_csv(cluster_markers_csv)

    col_gene = "names" if "names" in cm.columns else ("gene" if "gene" in cm.columns else None)
    col_group = "group" if "group" in cm.columns else ("cluster" if "cluster" in cm.columns else None)

    if col_gene is None or col_group is None:
        logging.warning("Could not find expected columns ('group'/'names') in cluster_markers CSV; skipping.")
        return

    cm[col_group] = cm[col_group].astype(str)
    subset = cm[(cm[col_group] == str(target_cluster)) & (cm[col_gene].isin(tls_genes))].copy()

    if subset.empty:
        logging.warning(f"No TLS genes found in DE list for cluster {target_cluster}.")
        return

    if "pvals_adj" in subset.columns:
        subset = subset.sort_values("pvals_adj")
    elif "scores" in subset.columns:
        subset = subset.sort_values("scores", ascending=False)

    out_path = os.path.join(outdir, f"cluster_{target_cluster}_TLS_gene_DE_check.csv")
    subset.to_csv(out_path, index=False)
    logging.info(f"Saved TLS gene DE check: {out_path}")


def find_spatially_variable_genes_moran(adata, outdir, n_hvgs=2000, top_n=200):
    """
    Optional: spatially variable genes (Moran's I) over HVGs.
    Results stored in adata.uns["moranI"] in many Squidpy versions.
    """
    logging.info("Running HVG selection for SVG (Moran's I)...")
    if "highly_variable" not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_hvgs)

    hvgs = adata.var_names[adata.var["highly_variable"]].tolist()
    if not hvgs:
        logging.warning("No HVGs found; skipping SVG Moran.")
        return

    logging.info("Computing Moran's I for HVGs...")
    sq.gr.spatial_autocorr(adata, genes=hvgs, mode="moran", n_perms=100)

    moran_df = adata.uns["moranI"].copy().sort_values("I", ascending=False)
    if "gene" not in moran_df.columns:
        moran_df.insert(0, "gene", moran_df.index)

    num_cols = moran_df.select_dtypes(include=[np.number]).columns
    moran_df[num_cols] = moran_df[num_cols].round(3)

    out_all = os.path.join(outdir, "moranI_hvgs_all.csv")
    out_top = os.path.join(outdir, f"moranI_hvgs_top{top_n}.csv")
    moran_df.to_csv(out_all, index=False)
    moran_df.head(top_n).to_csv(out_top, index=False)

    logging.info(f"Saved SVG Moran results: {out_all}")
    logging.info(f"Saved top {top_n} SVGs: {out_top}")


# =============================================================================
# Main
# =============================================================================

def main(args):
    """Run TLS-like spatial validation and export plots/CSVs."""
    os.makedirs(args.outdir, exist_ok=True)
    _setup_scanpy(dpi=args.dpi)
    sc.settings.figdir = args.outdir  # scanpy save= outputs go here

    logging.info(f"Output dir: {args.outdir}")
    logging.info("Loading AnnData...")
    adata = sc.read_h5ad(args.adata_path)
    logging.info(f"Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

    if args.cluster_key not in adata.obs:
        raise ValueError(f"[ERROR] '{args.cluster_key}' not found in adata.obs")
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str)

    if "spatial" not in adata.obsm:
        raise ValueError("[ERROR] adata.obsm['spatial'] not found. Required for Visium spatial plotting/stats.")

    # Marker sets (human, Visium)
    MARKERS = {
        "B_cells": ["CD19", "MS4A1", "CD79A", "CD74", "HLA-DRA"],
        "Plasma_cells": ["IGKC", "JCHAIN", "MZB1", "IGHG1"],
        "T_cells": ["TRAC", "CD3D", "CD3E", "IL7R", "LTB"],
        "CD8_T": ["CD8A", "CD8B", "GZMA", "NKG7"],
        "NK": ["NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "KLRK1"],
        "TLS_chemokines": ["CCL19", "CCL21", "CXCL13", "LTB"],
        "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    }

    # Module scores (obs columns named exactly as module keys)
    MODULES = {
        "B": ["CD19", "MS4A1", "CD79A", "CD74"],
        "Plasma": ["IGKC", "JCHAIN", "MZB1"],
        "T": ["TRAC", "CD3D", "CD3E", "IL7R"],
        "TLS": ["CCL19", "CCL21", "CXCL13", "LTB"],
        "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    }

    # Keep only genes present in the dataset
    MARKERS = {k: present_genes(adata, v) for k, v in MARKERS.items()}
    MODULES = {k: present_genes(adata, v) for k, v in MODULES.items()}

    logging.info("Present marker counts:")
    for k, v in MARKERS.items():
        logging.info(f"  - {k:15s}: {len(v)} genes")

    logging.info("Present module counts:")
    for k, v in MODULES.items():
        logging.info(f"  - {k:15s}: {len(v)} genes")

    # Compute module scores
    logging.info("Computing module scores...")
    for module_name, gene_list in MODULES.items():
        if len(gene_list) < 2:
            logging.warning(f"Skipping module '{module_name}' (needs >=2 present genes): {gene_list}")
            continue
        sc.tl.score_genes(adata, gene_list=gene_list, score_name=module_name, use_raw=False)

    score_cols = [m for m in MODULES.keys() if m in adata.obs.columns]
    logging.info(f"Scores computed: {score_cols}")

    # Spatial neighbors graph (needed for Moran + nhood enrichment)
    logging.info("Building spatial neighbors graph...")
    sq.gr.spatial_neighbors(adata, coord_type="grid")  # Visium
    logging.info("Spatial neighbors computed")

    # Optional: DE check for TLS genes in cluster markers table
    tls_check_genes = present_genes(
        adata,
        list(dict.fromkeys(MARKERS["TLS_chemokines"] + MARKERS["B_cells"] + MARKERS["T_cells"] + MARKERS["Plasma_cells"]))
    )
    export_tls_de_check(args.cluster_markers_csv, args.outdir, args.target_cluster, tls_check_genes)

    # Diagnostics: nonzero counts
    logging.info("Nonzero counts for TLS chemokines:")
    for g in ["CCL19", "CXCL13", "CCL21", "LTB"]:
        nz = count_nonzero(adata, g)
        logging.info(f"  - {g}: {'not in var_names' if nz is None else nz}")

    # Cluster localization plots
    logging.info("Saving cluster localization plots...")
    spatial_tissue_vs_clusters(
        adata,
        cluster_key=args.cluster_key,
        outdir=args.outdir,
        fname="spatial_clusters_overview.png",
        dots_marker_size=args.dots_marker_size,
        dpi=args.dpi,
    )
    spatial_tissue_vs_clusters(
        adata,
        cluster_key=args.cluster_key,
        outdir=args.outdir,
        fname=f"spatial_cluster_{args.target_cluster}_localization.png",
        target_cluster=args.target_cluster,
        dots_marker_size=args.dots_marker_size,
        dpi=args.dpi,
    )

    # Gene localization plots
    logging.info("Saving gene localization plots...")
    genes_to_plot = []
    for panel in ["TLS_chemokines", "B_cells", "T_cells", "Plasma_cells", "Epithelial"]:
        genes_to_plot.extend(MARKERS.get(panel, []))
    genes_to_plot = list(dict.fromkeys(genes_to_plot))  # de-dup preserve order

    for g in genes_to_plot:
        spatial_tissue_vs_dots(
            adata,
            g,
            args.outdir,
            fname_prefix="spatial_gene_localization",
            dots_marker_size=args.dots_marker_size,
            vmin=args.vmin_gene,
            vmax=args.vmax_gene,
            dpi=args.dpi,
        )

    # Module enrichment plots
    logging.info("Saving module enrichment plots...")
    for module_name in score_cols:
        spatial_tissue_vs_dots(
            adata,
            module_name,
            args.outdir,
            fname_prefix="spatial_module_enrichment",
            dots_marker_size=args.dots_marker_size,
            vmin=args.vmin_score,
            vmax=args.vmax_score,
            dpi=args.dpi,
        )

    # Moran's I (targeted) for TLS-related genes
    logging.info("Computing Moran's I for TLS-related genes...")
    moran_genes = present_genes(
        adata,
        ["CCL19", "CXCL13", "CCL21", "LTB",
         "TRAC", "CD3D",
         "CD79A", "MS4A1",
         "IGKC", "JCHAIN", "MZB1"],
    )

    if moran_genes:
        sq.gr.spatial_autocorr(adata, genes=moran_genes, mode="moran", n_perms=200)
        moran_df = adata.uns["moranI"].copy()
        if "gene" not in moran_df.columns:
            moran_df.insert(0, "gene", moran_df.index)

        num_cols = moran_df.select_dtypes(include=[np.number]).columns
        moran_df[num_cols] = moran_df[num_cols].round(3)

        moran_out = os.path.join(args.outdir, "moranI_tls_genes.csv")
        moran_df.to_csv(moran_out, index=False)
        logging.info(f"Saved Moran's I results: {moran_out}")
    else:
        logging.warning("No TLS genes present for Moran's I computation.")

    # Neighborhood enrichment (cluster adjacency)
    logging.info("Computing neighborhood enrichment (cluster adjacency)...")
    sq.gr.nhood_enrichment(
        adata,
        cluster_key=args.cluster_key,
        n_perms=args.nhood_n_perms,
        n_jobs=args.nhood_n_jobs,
    )

    # NOTE: squidpy save= writes into sc.settings.figdir and prefixes with leading "_"
    sq.pl.nhood_enrichment(
        adata,
        cluster_key=args.cluster_key,
        show=False,
        save="_nhood_enrichment.png",
    )
    logging.info("Saved neighborhood enrichment plot (_nhood_enrichment.png)")

    # Export z-scores if available (version-dependent)
    try:
        z = adata.uns["nhood_enrichment"]["zscore"]
        cats = adata.obs[args.cluster_key].astype("category").cat.categories
        z_df = pd.DataFrame(z, index=cats, columns=cats).round(3)
        z_out = os.path.join(args.outdir, "nhood_enrichment_zscores.csv")
        z_df.to_csv(z_out)
        logging.info(f"Saved: {z_out}")
    except Exception:
        logging.info("Could not export neighborhood enrichment z-scores (version-dependent). Plot saved.")

    # Export target cluster vs rest module score summary
    logging.info("Exporting target cluster vs rest module score summary...")
    mask_target = (adata.obs[args.cluster_key] == str(args.target_cluster))

    rows = []
    for module_name in score_cols:
        s0 = adata.obs.loc[mask_target, module_name].values
        sr = adata.obs.loc[~mask_target, module_name].values
        rows.append({
            "module": module_name,
            "target_mean": float(np.mean(s0)),
            "target_median": float(np.median(s0)),
            "rest_mean": float(np.mean(sr)),
            "rest_median": float(np.median(sr)),
            "delta_mean": float(np.mean(s0) - np.mean(sr)),
        })

    summary_df = pd.DataFrame(rows)
    num_cols = summary_df.select_dtypes(include=[np.number]).columns
    summary_df[num_cols] = summary_df[num_cols].round(3)
    summary_df = summary_df.sort_values("delta_mean", ascending=False)

    summary_out = os.path.join(args.outdir, f"cluster_{args.target_cluster}_score_summary.csv")
    summary_df.to_csv(summary_out, index=False)
    logging.info(f"Saved: {summary_out}")

    # Export spot-level module scores
    logging.info("Saving spot-level module scores...")
    spot_scores = adata.obs[[args.cluster_key] + score_cols].copy()
    spot_scores.insert(0, "spot_id", adata.obs_names)

    num_cols = spot_scores.select_dtypes(include=[np.number]).columns
    spot_scores[num_cols] = spot_scores[num_cols].round(3)

    spot_scores = spot_scores.sort_values(args.cluster_key)
    out_spot_scores = os.path.join(args.outdir, "spot_module_scores.csv")
    spot_scores.to_csv(out_spot_scores, index=False)
    logging.info(f"Saved: {out_spot_scores}")

    # Optional genome-wide SVG Moran on HVGs
    if args.run_svg:
        logging.info("Running optional HVG-wide Moran's I (SVG discovery)...")
        find_spatially_variable_genes_moran(adata, args.outdir, n_hvgs=args.hvg_n, top_n=args.svg_top_n)

    logging.info("DONE ✅ Outputs written to:")
    logging.info(f"  {args.outdir}/")
    logging.info("Key outputs:")
    logging.info("  - spatial_clusters_overview.png")
    logging.info(f"  - spatial_cluster_{args.target_cluster}_localization.png")
    logging.info("  - spatial_gene_localization_*.png")
    logging.info("  - spatial_module_enrichment_*.png")
    logging.info("  - moranI_tls_genes.csv")
    logging.info(f"  - cluster_{args.target_cluster}_score_summary.csv")
    logging.info("  - spot_module_scores.csv")
    logging.info("  - _nhood_enrichment.png (in outdir via scanpy figdir)")


if __name__ == "__main__":
    setup_logging()
    args = parse_args()
    main(args)
