#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import sys
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from umap.umap_ import UMAP
import matplotlib.pyplot as plt


INPUT_CSV      = "admet_filtered.csv"
OUTPUT_ALL     = "admet_umap_clusters.csv"
OUTPUT_MEDOIDS = "cluster_medoids.csv"
RAW_PLOT       = "umap_raw.png"
CLUSTER_PLOT   = "umap_clusters_medoids.png"


FEATURES = [
    'logS','logP','TPSA','nRot','BBB','MW','nHD','nHA',
    'logVDss','cl-plasma','t0.5','PPB','Fsp3'
]
K = 12


def main():
    if not os.path.exists(INPUT_CSV):
        sys.exit(f"Error: {INPUT_CSV} not found")
    df = pd.read_csv(INPUT_CSV)


    scaler = StandardScaler()
    X = scaler.fit_transform(df[FEATURES])


    umap_model = UMAP(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
    coords = umap_model.fit_transform(X)
    df['UMAP1'], df['UMAP2'] = coords[:,0], coords[:,1]
    df.to_csv(OUTPUT_ALL, index=False)


    # raw plot (same aspect ratio as cluster)
    plt.figure(figsize=(8,5))
    plt.scatter(df['UMAP1'], df['UMAP2'], color='gray', s=10, alpha=0.6)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(RAW_PLOT, dpi=300)
    plt.close()


    km = KMeans(n_clusters=K, random_state=42)
    df['cluster'] = km.fit_predict(coords)


    medoids = []
    for cid in range(K):
        subset = df[df['cluster'] == cid]
        center = km.cluster_centers_[cid]
        dists = np.linalg.norm(subset[['UMAP1','UMAP2']].values - center, axis=1)
        medoids.append(subset.iloc[dists.argmin()])
    medoids_df = pd.DataFrame(medoids)
    medoids_df.to_csv(OUTPUT_MEDOIDS, index=False)


    # print medoids
    for i, row in medoids_df.reset_index(drop=True).iterrows():
        print(f"Medoid {i+1}: {row['smiles']}")


    plt.figure(figsize=(10,5))
    scatter = plt.scatter(
        df['UMAP1'], df['UMAP2'],
        c=df['cluster'], cmap='tab20',
        s=30, alpha=0.8
    )
    plt.scatter(
        medoids_df['UMAP1'], medoids_df['UMAP2'],
        facecolors='none', edgecolors='black',
        s=120, marker='X', linewidth=1.5
    )
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")


    handles, _ = scatter.legend_elements(prop="colors", num=K)
    labels = [str(i+1) for i in range(K)]
    medoid_handle = plt.Line2D(
        [0], [0], marker='X', color='black',
        linestyle='None', markersize=10, markeredgewidth=1.5
    )
    handles.append(medoid_handle)
    labels.append('Medoids')


    plt.legend(
        handles, labels,
        title="Cluster",
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        fontsize=8,
        title_fontsize=10
    )
    plt.tight_layout(rect=(0,0,0.8,1))
    plt.savefig(CLUSTER_PLOT, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()



