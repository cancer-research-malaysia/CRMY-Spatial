import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sc
import stereo as st
import warnings
warnings.filterwarnings('ignore')
from scipy.io import mmwrite

datalist = [
    "SD507","SD560","SD1043", # A3B del/del
    "SD1182","SD1225", # A3B del/WT
    "SD683","SD693","SD781" # WT
]
bin_size = 30
n_pcs = 20
n_neighbors = 10

datas = []
for d in datalist:
    if d in ["SD507","SD693"]:
        data_path = "/stereoseq/all_samples/gef/SS200000873BR_D5.tissue.gef" # 100
    elif d in ["SD560","SD683"]:
        data_path = "/stereoseq/all_samples/gef/B03601F5.tissue.gef" # 100
    elif d in ["SD1043","SD1225"]:
        data_path = "/stereoseq/all_samples/gef/D04054A3.tissue.gef" # 100
    elif d in ["SD781","SD1182"]:
        data_path = "/stereoseq/all_samples/gef/D04054A6.tissue.gef" # 100
        
    data = st.io.read_gef(file_path=data_path, bin_size=bin_size)
    fil = np.array(pd.read_csv("filter1.csv")["gene"]) # MT, RPL, RPS
    data = data.tl.filter_genes(gene_list=fil)
    # sub by index
    #ids = pd.read_excel("bin{}_index.xlsx".format(bin_size), d, engine="openpyxl")
    ids = pd.read_csv("/stereoseq/all_samples/misc/{}.bin{}.coords.csv".format(d, bin_size))
    data = data.sub_by_index(ids.sort_values(by="bin")["index"])
    data.tl.cal_qc()
    data.tl.filter_cells(min_counts=20, min_genes=3, pct_counts_mt=5, inplace=True)
    data.tl.raw_checkpoint()
    data.tl.normalize_total(target_sum=10000) # normalise target
    data.tl.log1p() # transform exp matrix log
    data.tl.highly_variable_genes( # highly variable genes
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            n_top_genes=2000,
            res_key='highly_variable_genes'
            )
    data.tl.scale(max_value=10, zero_center=True) # scale"""
    """data.tl.pca( # pca
            use_highly_genes=False,
            n_pcs=30,
            res_key='pca'
            )
    data.tl.neighbors( # neighbours
            pca_res_key='pca',
            n_pcs=30,
            res_key='neighbors'
            )
    data.tl.umap( # umap
            pca_res_key='pca', 
            neighbors_res_key='neighbors', 
            res_key='umap',
            min_dist=0.5, spread=5
            )
    data.tl.leiden( # leiden clustering neighbours
            neighbors_res_key='neighbors', 
            res_key='leiden'
            )
    data.tl.louvain( # louvain clustering neighbours
            neighbors_res_key='neighbors', 
            res_key='louvain'
            )
    data.tl.phenograph( # phenograph clustering
            phenograph_k=30, 
            pca_res_key='pca', 
            res_key='phenograph'
            )
    data.tl.find_marker_genes( # find marker genes
            cluster_res_key='leiden',
            method='t_test',
            use_highly_genes=False,
            use_raw=True,
            res_key='t_marker_genes'
            )
    data.tl.filter_marker_genes( # filter marker genes
        marker_genes_res_key='t_marker_genes',
        min_fold_change=1,
        min_in_group_fraction=0.25,
        max_out_group_fraction=0.5,
        res_key='t_marker_genes_filtered'
    )"""

    #df = data.to_df()

    # export data into h5ad
    #data_path = sample_path + "h5ad/B03601F5.bin100.{}.processed.h5ad".format(sample)
    #data_path = "/stereoseq/all_samples/h5ad/bin{0}/{1}.bin{0}.processed.h5ad".format(bin_size,sample)
    #data_path = "/stereoseq/all_samples/h5ad/bin{0}/{1}.bin{0}.h5ad".format(bin_size,d)
    #st.io.write_h5ad(data, output=data_path)
    out_path = "/stereoseq/all_samples/h5ad/bin{0}/{1}.bin{0}.seurat.h5ad".format(bin_size,d)
    #adata = st.io.stereo_to_anndata(data, output=out_path, flavor="seurat")

    """df1 = pd.DataFrame(index=df.index)
    df1["x"] = data.position[:,0]
    df1["y"] = data.position[:,1]
    df1.reset_index(inplace=True)"""
    #df1.to_csv("/stereoseq/all_samples/h5ad/bin{0}/{1}.bin{0}.positions.csv".format(bin_size,d))

    #df.index.name = "barcodes"
    #df.to_csv("/stereoseq/all_samples/h5ad/bin{0}/{1}.bin{0}.processed.counts.csv".format(bin_size,d))

    target = "/stereoseq/all_samples/h5ad/bin{0}/{1}.bin{0}.raw.counts.int.mtx".format(bin_size,d)
    mmwrite(target, data.raw._exp_matrix.astype(int))
