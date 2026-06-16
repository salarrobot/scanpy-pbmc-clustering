from __future__ import annotations

import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import anndata
import scanpy as sc

FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(FIGDIR, exist_ok=True)
sc.settings.verbosity = 1
sc.settings.n_jobs = 1
sc.settings.figdir = FIGDIR
sc.set_figure_params(dpi=120, facecolor="white")
np.random.seed(0)

CELL_TYPES = ["CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells",
              "NK cells", "FCGR3A+ Monocytes", "Dendritic cells", "Megakaryocytes"]



def problem1():
    print("\n" + "#" * 70)
    print("# PROBLEM 1 - Clustering 3k PBMCs")
    print("#" * 70)

    adata = sc.datasets.pbmc3k()             
    adata.var_names_make_unique()

    # Basic filtering + QC
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    adata = adata[(adata.obs.n_genes_by_counts < 2500)
                  & (adata.obs.n_genes_by_counts > 200)
                  & (adata.obs.pct_counts_mt < 5), :].copy()
    adata.layers["counts"] = adata.X.copy()
    print("After QC:", adata.shape)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, layer="counts", n_top_genes=2000,
                                min_mean=0.0125, max_mean=3, min_disp=0.5,
                                flavor="seurat_v3")

    adata.layers["scaled"] = adata.X.toarray()
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"], layer="scaled")
    sc.pp.scale(adata, max_value=10, layer="scaled")

    # PCA / neighbors / UMAP / Leiden
    sc.pp.pca(adata, layer="scaled", svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.7, random_state=0,
                 flavor="igraph", n_iterations=2, directed=False)

    n_clusters = adata.obs["leiden"].cat.categories.size
    print("Number of Leiden clusters:", n_clusters)

    markers_a = ["NKG7", "IL7R", "CD14", "MS4A1", "CD8A", "FCGR3A", "FCER1A", "PPBP"]
    sc.pl.umap(adata, color=["leiden"] + markers_a, ncols=3, frameon=False,
               cmap="viridis", show=False, save="_p1a_markers.png")
    print("(a) Saved UMAP of Leiden clusters + markers -> figures/umap_p1a_markers.png")

    verify = {"CD4 T cells": ["IL7R"], "CD14+ Monocytes": ["CD14", "LYZ"],
              "B cells": ["MS4A1", "CD79A"], "CD8 T cells": ["CD8A"],
              "NK cells": ["GNLY", "NKG7"], "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7"],
              "Dendritic cells": ["FCER1A", "CST3"], "Megakaryocytes": ["PPBP"]}
    genes = [g for gs in verify.values() for g in gs if g in adata.var_names]
    mean_expr = sc.get.obs_df(adata, keys=genes + ["leiden"]
                              ).groupby("leiden", observed=True).mean()

    tutorial = ["CD4 T cells", "B cells", "CD14+ Monocytes", "NK cells",
                "CD8 T cells", "FCGR3A+ Monocytes", "Dendritic cells", "Megakaryocytes"]
    cats = list(adata.obs["leiden"].cat.categories)
    if n_clusters == len(tutorial):
        mapping = dict(zip(cats, tutorial))
    else:
        score = pd.DataFrame(index=mean_expr.index)
        for ct, gs in verify.items():
            gs = [g for g in gs if g in mean_expr.columns]
            z = (mean_expr[gs] - mean_expr[gs].mean()) / (mean_expr[gs].std() + 1e-9)
            score[ct] = z.mean(axis=1)
        mapping = score.idxmax(axis=1).to_dict()

    adata.obs["cell_type"] = adata.obs["leiden"].map(mapping).astype("category")
    sc.pl.umap(adata, color="cell_type", legend_loc="on data", frameon=False,
               title="Cell types", show=False, save="_p1b_celltypes.png")

    counts = adata.obs["cell_type"].value_counts().reindex(CELL_TYPES).astype(int)
    print("\n(b) Number of cells per cell type:")
    print(counts.to_string())
    print("Total:", int(counts.sum()))



def problem2():
    print("\n" + "#" * 70)
    print("# PROBLEM 2 - Data integration (ingest & BBKNN)")
    print("#" * 70)

    adata_ref = sc.datasets.pbmc3k_processed()      
    adata = sc.datasets.pbmc68k_reduced()           
    var_names = adata_ref.var_names.intersection(adata.var_names)
    adata_ref = adata_ref[:, var_names].copy()
    adata = adata[:, var_names].copy()
    print("Reference:", adata_ref.shape, "| Query:", adata.shape,
          "| shared genes:", len(var_names))

    ref_b, query_b = adata_ref.copy(), adata.copy()   

    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)
    sc.tl.ingest(adata, adata_ref, obs="louvain")
    adata.uns["louvain_colors"] = adata_ref.uns["louvain_colors"]
    sc.pl.umap(adata, color=["louvain", "bulk_labels"], wspace=0.5,
               show=False, save="_p2a_query_ingest.png")

    ingest_counts = (adata.obs["louvain"].value_counts()
                     .reindex(CELL_TYPES).fillna(0).astype(int))
    print("\n(a) ingest - number of query cells per cell type:")
    print(ingest_counts.to_string())
    print("Total:", int(ingest_counts.sum()))

    adata_concat = anndata.concat([ref_b, query_b], label="batch", keys=["ref", "new"])
    adata_concat.obs["batch"] = adata_concat.obs["batch"].astype("category")
    sc.tl.pca(adata_concat)
    sc.external.pp.bbknn(adata_concat, batch_key="batch")
    sc.tl.umap(adata_concat)
    sc.pl.umap(adata_concat, color=["batch"], show=False,
               save="_p2b_concat_bbknn_batch.png")

    
    cats = list(ref_b.obs["louvain"].cat.categories)
    ref_mask = (adata_concat.obs["batch"].values == "ref")
    onehot = np.zeros((adata_concat.n_obs, len(cats)))
    onehot[np.where(ref_mask)[0], ref_b.obs["louvain"].cat.codes.values] = 1.0
    votes = adata_concat.obsp["connectivities"].dot(onehot)
    q = np.where(~ref_mask)[0]
    no_vote = votes[q].sum(axis=1) == 0
    if no_vote.any():                                
        votes2 = adata_concat.obsp["connectivities"].dot(votes)
        votes[q[no_vote]] = votes2[q[no_vote]]
    pred = np.array(cats)[votes[q].argmax(axis=1)]
    query_b.obs["louvain_bbknn"] = pd.Categorical(pred, categories=cats)

    bbknn_counts = (query_b.obs["louvain_bbknn"].value_counts()
                    .reindex(CELL_TYPES).fillna(0).astype(int))
    print("\n(b) BBKNN - number of query cells per cell type:")
    print(bbknn_counts.to_string())
    print("Total:", int(bbknn_counts.sum()))

    print("\nComparison (700 query cells):")
    print(pd.DataFrame({"ingest (a)": ingest_counts,
                        "bbknn (b)": bbknn_counts}).to_string())


if __name__ == "__main__":
    problem1()
    problem2()
    print("\nDone. Figures saved in:", FIGDIR)
