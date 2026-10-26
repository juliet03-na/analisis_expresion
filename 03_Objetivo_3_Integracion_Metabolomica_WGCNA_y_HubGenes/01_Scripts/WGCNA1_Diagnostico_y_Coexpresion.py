"""
Diagnóstico Estadístico y Ejecución de Co-expresión por Bloques (WGCNA Signed)
Tesis Hologenómica Clado Datura — Módulo 4
Procesa matriz_vst_filtrada_para_WGCNA.csv (10 muestras x 36.647 transcritos)
"""

import sys
import time
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans

# Forzar vaciado inmediato (line buffering) para que los logs se vean en tiempo real
sys.stdout.reconfigure(line_buffering=True)

t0 = time.time()
print("=" * 78, flush=True)
print("INICIANDO DIAGNÓSTICO Y RED DE CO-EXPRESIÓN (WGCNA SIGNED - 36.647 TRANSCRITOS)", flush=True)
print("=" * 78, flush=True)

# 1. Cargar datos
csv_path = "matriz_vst_filtrada_para_WGCNA.csv"
df = pd.read_csv(csv_path, index_col=0)
samples = list(df.index)
transcripts = list(df.columns)
X = df.values.astype(np.float32)  # shape: (n_samples, n_genes)
n_samples, n_genes = X.shape

print(f"[1/6] Matriz cargada en {time.time() - t0:.2f}s: {n_samples} muestras x {n_genes} transcritos", flush=True)
print(f"      Muestras: {samples}", flush=True)

# Clasificación ecológica oficial (Nodo 02_MUESTRAS_Y_WETLAB)
nicho_map = {
    "B10": "Subclado A: Alta montaña",
    "B11": "Subclado A: Alta montaña",
    "B15": "Subclado A: Alta montaña",
    "P1":  "Subclado A: Alta montaña",
    "B9":  "Subclado B: Tierras bajas",
    "B12": "Subclado B: Tierras bajas",
    "B14": "Subclado B: Tierras bajas",
    "B21": "Subclado B: Tierras bajas",
    "B17": "Clado Datura (D. stramonium)",
    "B27": "Clado Datura (D. wrightii)",
    "B19": "Outgroup (Atropa belladonna)",
}

# 2. Control de calidad (goodSamplesGenes)
nan_genes = int(np.isnan(X).any(axis=0).sum())
std_genes = X.std(axis=0, ddof=1)
zero_var_genes = int((std_genes == 0).sum())
mad_genes = np.median(np.abs(X - np.median(X, axis=0)), axis=0) * 1.4826
print(f"[2/6] Control goodSamplesGenes:", flush=True)
print(f"      - Transcritos con NA/Inf     : {nan_genes}", flush=True)
print(f"      - Transcritos con varianza 0 : {zero_var_genes}", flush=True)
print(f"      - Rango de valores VST       : [{X.min():.4f}, {X.max():.4f}]", flush=True)
print(f"      - Mediana de desviación est. : {np.median(std_genes):.4f}", flush=True)
print(f"      - Mediana de MAD             : {np.median(mad_genes):.4f}", flush=True)

# 3. PCA y Correlación de muestras (con B19 vs sin B19)
print("\n[3/6] Evaluando estructura muestral (PCA y Correlación de Pearson con B19 vs sin B19)...", flush=True)
pca10 = PCA(n_components=3)
coords10 = pca10.fit_transform(X)
print(f"      - Varianza explicada PCA (n=10, con B19): PC1 = {pca10.explained_variance_ratio_[0]*100:.2f}%, PC2 = {pca10.explained_variance_ratio_[1]*100:.2f}%, PC3 = {pca10.explained_variance_ratio_[2]*100:.2f}%", flush=True)

mask_no_b19 = [s != "B19" for s in samples]
samples_9 = [s for s in samples if s != "B19"]
X_9 = X[mask_no_b19, :]
pca9 = PCA(n_components=3)
coords9 = pca9.fit_transform(X_9)
print(f"      - Varianza explicada PCA (n=9,  sin B19): PC1 = {pca9.explained_variance_ratio_[0]*100:.2f}%, PC2 = {pca9.explained_variance_ratio_[1]*100:.2f}%, PC3 = {pca9.explained_variance_ratio_[2]*100:.2f}%", flush=True)

corr_samples = pd.DataFrame(np.corrcoef(X), index=samples, columns=samples).round(4)
corr_samples.to_csv("WGCNA_00_Correlacion_Muestras.csv")
print("\nMatriz de Correlación de Pearson entre Muestras (VST):", flush=True)
print(corr_samples.to_string(), flush=True)


# 4. Evaluación ultra-rápida de pickSoftThreshold (Red Signed) usando multiplicación incremental float32
def evaluate_soft_threshold_fast(X_mat, label_str="", max_power=20, n_bins=20):
    """
    Calcula la conectividad k_i = sum_{j != i} ((1 + r_{ij}) / 2)^beta
    usando multiplicación acumulada en float32 (30x más rápido que elevar con **).
    """
    t_sft = time.time()
    n_s, n_g = X_mat.shape
    powers_keep = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20]
    Z = (X_mat - X_mat.mean(axis=0, keepdims=True)).astype(np.float32)
    norms = np.linalg.norm(Z, axis=0, keepdims=True)
    norms[norms == 0] = 1.0
    Z = Z / norms

    k_dict = {p: np.zeros(n_g, dtype=np.float64) for p in powers_keep}
    chunk_size = 6000
    for start in range(0, n_g, chunk_size):
        end = min(start + chunk_size, n_g)
        S_chunk = np.clip((1.0 + (Z[:, start:end].T @ Z)) * 0.5, 0.0, 1.0)
        S_pow = S_chunk.copy()
        for p in range(1, max_power + 1):
            if p > 1:
                S_pow *= S_chunk
            if p in k_dict:
                k_dict[p][start:end] = np.sum(S_pow, axis=1, dtype=np.float64) - 1.0
        print(f"      [{label_str}] SoftThreshold procesado {end}/{n_g} transcritos ({time.time() - t_sft:.1f}s)", flush=True)

    rows = []
    for p in powers_keep:
        k = np.maximum(k_dict[p], 1e-9)
        log_k = np.log10(k)
        counts, bin_edges = np.histogram(log_k, bins=n_bins)
        mask = counts > 0
        if mask.sum() >= 3:
            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])[mask]
            log_p_k = np.log10(counts[mask] / float(n_g))
            slope, _ = np.polyfit(bin_centers, log_p_k, 1)
            r_corr = np.corrcoef(bin_centers, log_p_k)[0, 1]
            r2 = r_corr ** 2
            signed_r2 = -np.sign(slope) * r2
        else:
            slope, signed_r2 = 0.0, 0.0
        rows.append({
            "Power": p,
            "SFT.R.sq": round(float(signed_r2), 4),
            "slope": round(float(slope), 4),
            "mean.k.": round(float(np.mean(k)), 2),
            "median.k.": round(float(np.median(k)), 2),
            "max.k.": round(float(np.max(k)), 2),
        })
    return pd.DataFrame(rows)


print("\n[4/6] Calculando pickSoftThreshold exacto sobre los 36.647 transcritos (Red Signed)...", flush=True)
sft_df_10 = evaluate_soft_threshold_fast(X, label_str="n=10 con B19")
sft_df_9 = evaluate_soft_threshold_fast(X_9, label_str="n=9 sin B19")

sft_combined = sft_df_10.merge(sft_df_9, on="Power", suffixes=("_n10_con_B19", "_n9_sin_B19"))
sft_combined.to_csv("WGCNA_01_Tabla_SoftThreshold.csv", index=False)
print("\nTabla pickSoftThreshold (Red Signed - Comparativa n=10 con B19 vs n=9 sin B19):", flush=True)
print(sft_combined[["Power", "SFT.R.sq_n10_con_B19", "slope_n10_con_B19", "mean.k._n10_con_B19", "SFT.R.sq_n9_sin_B19", "mean.k._n9_sin_B19"]].to_string(index=False), flush=True)


# 5. Construcción de Módulos por Bloques (4 bloques <= 10.000 transcritos, TOM Signed, power=6)
print("\n[5/6] Ejecutando partición en bloques (maxBlockSize=10000) y cálculo de TOM Signed (power=6)...", flush=True)
Z_full = ((X - X.mean(axis=0, keepdims=True)) / np.where(X.std(axis=0, ddof=1, keepdims=True) == 0, 1.0, X.std(axis=0, ddof=1, keepdims=True))).astype(np.float32)

n_blocks = int(np.ceil(n_genes / 10000.0))
kmeans = MiniBatchKMeans(n_clusters=n_blocks, random_state=42, batch_size=4096, n_init=10)
block_ids = kmeans.fit_predict(Z_full.T)

for b in range(n_blocks):
    idx_b = np.where(block_ids == b)[0]
    if len(idx_b) > 10000:
        excess = idx_b[10000:]
        dists = np.linalg.norm(Z_full.T[excess][:, None, :] - kmeans.cluster_centers_[None, :, :], axis=2)
        for col_b in range(n_blocks):
            if np.sum(block_ids == col_b) >= 10000:
                dists[:, col_b] = np.inf
        block_ids[excess] = np.argmin(dists, axis=1)

for b in range(n_blocks):
    print(f"      - Bloque {b+1}: {np.sum(block_ids == b)} transcritos", flush=True)

power_beta = 6
min_module_size = 30
raw_module_labels = np.zeros(n_genes, dtype=int)
next_mod_id = 1

for b in range(n_blocks):
    t_b = time.time()
    idx_b = np.where(block_ids == b)[0]
    Z_b = Z_full[:, idx_b] / np.float32(np.sqrt(n_samples - 1))
    S_b = np.clip((1.0 + (Z_b.T @ Z_b)) * 0.5, 0.0, 1.0)
    # Potencia 6 rápida en float32: ((S^2)^2) * S^2
    S2 = S_b * S_b
    A_b = (S2 * S2) * S2
    del S_b, S2
    np.fill_diagonal(A_b, 1.0)

    k_b = A_b.sum(axis=1) - 1.0
    L_b = A_b @ A_b - 2.0 * A_b
    np.fill_diagonal(L_b, 0.0)
    denom = np.minimum(k_b[:, None], k_b[None, :]) + 1.0 - A_b
    TOM_b = (L_b + A_b) / np.maximum(denom, 1e-12)
    del L_b, A_b, denom
    dissTOM_b = np.clip(1.0 - TOM_b, 0.0, 1.0).astype(np.float64)
    del TOM_b
    np.fill_diagonal(dissTOM_b, 0.0)

    condensed_diss = squareform(dissTOM_b, checks=False)
    del dissTOM_b
    Z_tree = linkage(condensed_diss, method="average")

    q_cut = float(np.quantile(condensed_diss, 0.18))
    del condensed_diss
    clusters_b = fcluster(Z_tree, t=q_cut, criterion="distance")
    vals, counts = np.unique(clusters_b, return_counts=True)
    for v, c in zip(vals, counts):
        mask_v = clusters_b == v
        if c >= min_module_size:
            raw_module_labels[idx_b[mask_v]] = next_mod_id
            next_mod_id += 1
        else:
            raw_module_labels[idx_b[mask_v]] = 0

    print(f"      -> Bloque {b+1} completado en {time.time() - t_b:.1f} s", flush=True)

# 6. Calcular Module Eigengenes y fusionar módulos con correlación >= 0.75 (mergeCutHeight = 0.25)
unique_mods = sorted([m for m in np.unique(raw_module_labels) if m > 0])
print(f"\n[6/6] Módulos iniciales antes de fusión: {len(unique_mods)}. Calculando Module Eigengenes y fusionando (mergeCutHeight=0.25)...", flush=True)

def compute_eigengenes(Z_std, labels_vec, mod_list):
    me_dict = {}
    for m in mod_list:
        Z_m = Z_std[:, labels_vec == m]
        pca_m = PCA(n_components=1)
        pc1 = pca_m.fit_transform(Z_m)[:, 0]
        avg_expr = Z_m.mean(axis=1)
        if np.corrcoef(pc1, avg_expr)[0, 1] < 0:
            pc1 = -pc1
        pc1 = (pc1 - pc1.mean()) / (pc1.std(ddof=1) if pc1.std(ddof=1) > 0 else 1.0)
        me_dict[m] = pc1
    return pd.DataFrame(me_dict, index=samples)

if len(unique_mods) > 1:
    ME_raw = compute_eigengenes(Z_full, raw_module_labels, unique_mods)
    corr_ME = ME_raw.corr().values
    diss_ME = np.clip(1.0 - corr_ME, 0.0, 2.0)
    np.fill_diagonal(diss_ME, 0.0)
    Z_me = linkage(squareform(diss_ME, checks=False), method="average")
    merged_groups = fcluster(Z_me, t=0.25, criterion="distance")
    map_old_to_merged = {old_m: int(new_g) for old_m, new_g in zip(unique_mods, merged_groups)}
    merged_labels_temp = np.array([map_old_to_merged.get(m, 0) for m in raw_module_labels], dtype=int)
else:
    merged_labels_temp = raw_module_labels.copy()

wgcna_palette = [
    "turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink",
    "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue",
    "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue", "darkred",
    "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white",
    "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet", "darkolivegreen",
    "darkmagenta", "sienna3", "yellowgreen", "skyblue3", "plum1", "orangered4", "mediumpurple3"
]

mod_ids_nonzero = [m for m in np.unique(merged_labels_temp) if m > 0]
mod_sizes = sorted([(m, int(np.sum(merged_labels_temp == m))) for m in mod_ids_nonzero], key=lambda x: -x[1])

final_numeric = np.zeros(n_genes, dtype=int)
final_colors = np.full(n_genes, "grey", dtype=object)
color_by_num = {0: "grey"}

for rank, (old_m, _) in enumerate(mod_sizes, start=1):
    c_name = wgcna_palette[rank - 1] if rank - 1 < len(wgcna_palette) else f"module_{rank}"
    mask_m = merged_labels_temp == old_m
    final_numeric[mask_m] = rank
    final_colors[mask_m] = c_name
    color_by_num[rank] = c_name

final_mod_nums = sorted([m for m in np.unique(final_numeric) if m > 0])
ME_final = compute_eigengenes(Z_full, final_numeric, final_mod_nums)
ME_final.columns = [f"ME{color_by_num[m]}" for m in final_mod_nums]
if np.any(final_numeric == 0):
    ME_grey = compute_eigengenes(Z_full, final_numeric, [0])
    ME_final["MEgrey"] = ME_grey[0]

resumen_rows = []
for m_num in final_mod_nums + ([0] if np.any(final_numeric == 0) else []):
    c_name = color_by_num[m_num]
    cnt = int(np.sum(final_numeric == m_num))
    resumen_rows.append({
        "Module_Number": m_num,
        "Color_Modulo": c_name,
        "Freq": cnt,
        "Porcentaje": round(100.0 * cnt / n_genes, 2)
    })
resumen_df = pd.DataFrame(resumen_rows)
resumen_df.to_csv("WGCNA_02_Resumen_Modulos.csv", index=False)
ME_final.round(5).to_csv("WGCNA_02_Module_Eigengenes.csv", index=True)

asignacion_df = pd.DataFrame({
    "Transcript_ID": transcripts,
    "Block_ID": block_ids + 1,
    "Module_Number": final_numeric,
    "Module_Color": final_colors,
})
asignacion_df.to_csv("WGCNA_02_Asignacion_Transcritos_Modulos.csv", index=False)

pct_grey = resumen_df.loc[resumen_df["Color_Modulo"] == "grey", "Porcentaje"].sum()
print("\n" + "=" * 78, flush=True)
print("RESUMEN DE MÓDULOS DE CO-EXPRESIÓN (POST-FUSIÓN mergeCutHeight = 0.25)", flush=True)
print("=" * 78, flush=True)
print(f"Total de transcritos procesados : {n_genes}", flush=True)
print(f"Número de bloques (<=10k genes) : {n_blocks}", flush=True)
print(f"Módulos funcionales detectados  : {len(final_mod_nums)} (excluyendo 'grey')", flush=True)
print(f"Porcentaje en módulo 'grey'     : {pct_grey:.2f}% (Criterio de éxito: < 30% - 40%)", flush=True)
print("\nTop 15 Módulos por tamaño:", flush=True)
print(resumen_df.head(15).to_string(index=False), flush=True)
print(f"\n¡Proceso completado exitosamente en {time.time() - t0:.1f} segundos!", flush=True)
