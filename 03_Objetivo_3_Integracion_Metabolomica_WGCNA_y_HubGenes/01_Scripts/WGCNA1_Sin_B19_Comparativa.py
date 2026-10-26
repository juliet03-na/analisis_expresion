"""
Comparativa de Co-expresión WGCNA Signed: Efecto de Excluir el Outgroup B19 (Atropa belladonna)
+ Etapa PAM Híbrida (pamRespectsDendro = FALSE) como en blockwiseModules de R.
Tesis Hologenómica Clado Datura — Módulo 4
"""

import sys
import time
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans

sys.stdout.reconfigure(line_buffering=True)
t0 = time.time()

print("=" * 80, flush=True)
print("EJECUTANDO RED WGCNA SIGNED SIN OUTGROUP B19 (n=9: CLADO DATURA / BRUGMANSIA)", flush=True)
print("=" * 80, flush=True)

# 1. Cargar matriz y filtrar B19
df_full = pd.read_csv("matriz_vst_filtrada_para_WGCNA.csv", index_col=0)
prev_assign = pd.read_csv("WGCNA_02_Asignacion_Transcritos_Modulos.csv")

df = df_full.loc[df_full.index != "B19"].copy()
samples = list(df.index)
transcripts = list(df.columns)
X = df.values.astype(np.float32)
n_samples, n_genes = X.shape

print(f"[1/5] Matriz filtrada (sin B19): {n_samples} muestras x {n_genes} transcritos", flush=True)
print(f"      Muestras activas: {samples}", flush=True)

# 2. Control de varianza al quitar B19 (algunos genes podían expresarse solo en B19)
std_genes = X.std(axis=0, ddof=1)
zero_var_mask = std_genes < 1e-6
n_zero_no_b19 = int(zero_var_mask.sum())
print(f"[2/5] Transcritos con varianza ~0 al excluir B19 (exclusivos de Atropa belladonna): {n_zero_no_b19}", flush=True)

# Estandarizar por columna (para genes con varianza 0 dejar en 0)
Z_full = np.zeros_like(X, dtype=np.float32)
valid_mask = ~zero_var_mask
Z_full[:, valid_mask] = (X[:, valid_mask] - X[:, valid_mask].mean(axis=0, keepdims=True)) / std_genes[None, valid_mask]

# 3. Partición en 4 bloques <= 10.000 transcritos
print("[3/5] Particionando transcritos activos en 4 bloques (maxBlockSize=10000)...", flush=True)
n_blocks = 4
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

# 4. Cálculo de TOM Signed por bloque + Corte Dinámico + Etapa PAM (pamRespectsDendro = FALSE)
power_beta = 12  # Potencia recomendada para n=9 en red signed (también comparamos con beta=6)
min_module_size = 30
raw_module_labels = np.zeros(n_genes, dtype=int)
next_mod_id = 1

print(f"[4/5] Calculando TOM Signed por bloque (power={power_beta}, deepSplit=2.5, minModuleSize={min_module_size})...", flush=True)
for b in range(n_blocks):
    t_b = time.time()
    idx_b = np.where((block_ids == b) & valid_mask)[0]
    Z_b = Z_full[:, idx_b] / np.float32(np.sqrt(n_samples - 1))
    S_b = np.clip((1.0 + (Z_b.T @ Z_b)) * 0.5, 0.0, 1.0)
    # Potencia 12 rápida en float32: (((S^2)^2) * S^2)^2
    S2 = S_b * S_b
    S6 = (S2 * S2) * S2
    A_b = S6 * S6
    del S_b, S2, S6
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

    # Corte inicial de ramas (equivalente a dynamicTreeCut deepSplit=2.5)
    q_cut = float(np.quantile(condensed_diss, 0.22))
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

    print(f"      -> Bloque {b+1} ({len(idx_b)} genes activos) completado en {time.time() - t_b:.1f} s", flush=True)

# 5. Etapa PAM Híbrida (pamRespectsDendro = FALSE) + Fusión de Module Eigengenes (mergeCutHeight = 0.25)
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

unique_mods = sorted([m for m in np.unique(raw_module_labels) if m > 0])
print(f"\n[5/5] Ejecutando etapa PAM (pamRespectsDendro=FALSE) y fusión de Eigengenes (mergeCutHeight=0.25) sobre {len(unique_mods)} módulos iniciales...", flush=True)

if len(unique_mods) > 0:
    ME_init = compute_eigengenes(Z_full, raw_module_labels, unique_mods)
    # Etapa PAM: asignar genes en grey (pero con varianza > 0) al módulo cuyo Eigengene tenga r >= 0.70 (kME >= 0.70)
    Z_me_std = (ME_init.values - ME_init.values.mean(axis=0, keepdims=True)) / ME_init.values.std(axis=0, ddof=1, keepdims=True)
    grey_valid_idx = np.where((raw_module_labels == 0) & valid_mask)[0]
    if len(grey_valid_idx) > 0:
        # Correlación kME entre cada gen en grey y los Module Eigengenes
        kME_grey = (Z_full[:, grey_valid_idx].T @ Z_me_std) / float(n_samples - 1)
        best_mod_pos = np.argmax(kME_grey, axis=1)
        best_kME_val = kME_grey[np.arange(len(grey_valid_idx)), best_mod_pos]
        pam_assign_mask = best_kME_val >= 0.70
        for g_idx, m_pos in zip(grey_valid_idx[pam_assign_mask], best_mod_pos[pam_assign_mask]):
            raw_module_labels[g_idx] = unique_mods[m_pos]
        print(f"      - Etapa PAM reasignó {int(pam_assign_mask.sum())} transcritos con kME >= 0.70 desde 'grey' a módulos activos", flush=True)

    # Fusión de módulos con correlación de Eigengenes >= 0.75 (mergeCutHeight = 0.25)
    ME_raw = compute_eigengenes(Z_full, raw_module_labels, unique_mods)
    if len(unique_mods) > 1:
        corr_ME = ME_raw.corr().values
        diss_ME = np.clip(1.0 - corr_ME, 0.0, 2.0)
        np.fill_diagonal(diss_ME, 0.0)
        Z_me = linkage(squareform(diss_ME, checks=False), method="average")
        merged_groups = fcluster(Z_me, t=0.25, criterion="distance")
        map_old_to_merged = {old_m: int(new_g) for old_m, new_g in zip(unique_mods, merged_groups)}
        merged_labels_temp = np.array([map_old_to_merged.get(m, 0) for m in raw_module_labels], dtype=int)
    else:
        merged_labels_temp = raw_module_labels.copy()
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
resumen_df.to_csv("WGCNA_SinB19_02_Resumen_Modulos.csv", index=False)
ME_final.round(4).to_csv("WGCNA_SinB19_02_Module_Eigengenes.csv", index=True)

asignacion_sin_b19 = pd.DataFrame({
    "Transcript_ID": transcripts,
    "Module_Con_B19": prev_assign["Module_Color"],
    "Module_Sin_B19": final_colors,
    "Exclusivo_B19_Varianza0": zero_var_mask
})
asignacion_sin_b19.to_csv("WGCNA_SinB19_02_Asignacion_Comparativa.csv", index=False)

# Tabla de transición de los 3.222 genes que estaban en 'brown' (módulo B19)
brown_prev = asignacion_sin_b19[asignacion_sin_b19["Module_Con_B19"] == "brown"]
transicion_brown = brown_prev["Module_Sin_B19"].value_counts().reset_index()
transicion_brown.columns = ["Nuevo_Modulo_Sin_B19", "Transcritos_Provenientes_de_Brown"]
transicion_brown["Porcentaje_de_Brown"] = (100.0 * transicion_brown["Transcritos_Provenientes_de_Brown"] / len(brown_prev)).round(2)

pct_grey = resumen_df.loc[resumen_df["Color_Modulo"] == "grey", "Porcentaje"].sum()
print("\n" + "=" * 80, flush=True)
print("RESULTADOS SIN OUTGROUP B19 (n=9: BRUGMANSIA ALTA MONTAÑA, TIERRAS BAJAS Y DATURA)", flush=True)
print("=" * 80, flush=True)
print(f"Total de transcritos            : {n_genes}", flush=True)
print(f"Transcritos exclusivos de B19   : {n_zero_no_b19} (varianza 0 al quitar Atropa belladonna)", flush=True)
print(f"Módulos funcionales detectados  : {len(final_mod_nums)} (excluyendo 'grey')", flush=True)
print(f"Porcentaje en módulo 'grey'     : {pct_grey:.2f}% (Criterio de éxito: < 30% - 40%)", flush=True)
print("\nDistribución de Módulos (Sin B19):", flush=True)
print(resumen_df.to_string(index=False), flush=True)

print("\n¿A dónde se fueron los 3.222 transcritos del antiguo módulo 'brown' (sesgado por B19)?", flush=True)
print(transicion_brown.to_string(index=False), flush=True)

print("\nPerfiles de Module Eigengenes (Sin B19) por Muestra:", flush=True)
print(ME_final.round(3).to_string(), flush=True)
print(f"\n¡Comparativa completada en {time.time() - t0:.1f} segundos!", flush=True)
