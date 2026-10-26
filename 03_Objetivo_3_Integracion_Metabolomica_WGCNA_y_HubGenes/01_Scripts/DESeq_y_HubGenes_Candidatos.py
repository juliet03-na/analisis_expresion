"""
Integración Módulo 3 + Módulo 4:
1. Contraste de Expresión Diferencial: Subclado A (Alta Montaña: B10, B11, B15) vs Subclado B (Tierras Bajas: B9, B12, B14, B21)
2. Cálculo de Module Membership (kME) y Gene Significance (GS) contra los alcaloides clave (X70T, X52T, X82P, X71T)
3. Identificación de Hub Genes Candidatos (Intersección WGCNA + Expresión Diferencial, modelo híbrido Pablo/Ronald)
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)

# 1. Cargar matriz VST (36.647 transcritos), asignación de módulos (Sin B19), Eigengenes y Metaboloma homologado (n=8)
df_vst = pd.read_csv("matriz_vst_filtrada_para_WGCNA.csv", index_col=0)
df_mod = pd.read_csv("WGCNA_SinB19_02_Asignacion_Comparativa.csv")
df_me = pd.read_csv("WGCNA_SinB19_02_Module_Eigengenes.csv", index_col=0)
df_meta = pd.read_csv("fenotipos_TA_homologados_n8_SinOutgroups.csv", index_col=0)

alta_montana = ["B10", "B11", "B15"]
tierras_bajas = ["B9", "B12", "B14", "B21"]
datura_sp = ["B17", "B27"]

# 2. Expresión diferencial en escala log2 (VST): Alta Montaña (Subclado A) vs Tierras Bajas (Subclado B)
X_alta = df_vst.loc[alta_montana].values      # (3, 36647)
X_baja = df_vst.loc[tierras_bajas].values     # (4, 36647)
X_datura = df_vst.loc[datura_sp].values       # (2, 36647)

mean_alta = X_alta.mean(axis=0)
mean_baja = X_baja.mean(axis=0)
mean_datura = X_datura.mean(axis=0)

# En escala VST (que es asintóticamente log2), la diferencia de medias equivale al log2(Fold Change) encogido (shrunk LFC)
log2FC_Alta_vs_Baja = mean_alta - mean_baja
log2FC_Brugmansia_vs_Datura = df_vst.loc[alta_montana + tierras_bajas].values.mean(axis=0) - mean_datura

# Prueba t de Welch vectorizada entre Alta Montaña y Tierras Bajas
t_stat, p_val_de = stats.ttest_ind(X_alta, X_baja, axis=0, equal_var=False)
p_val_de = np.nan_to_num(p_val_de, nan=1.0)

# Ajuste FDR Benjamini-Hochberg
n_g = len(p_val_de)
order_p = np.argsort(p_val_de)
ranks = np.empty_like(order_p)
ranks[order_p] = np.arange(1, n_g + 1)
fdr_de = np.minimum(1.0, p_val_de * n_g / ranks)
# Asegurar monotonicidad en FDR
fdr_sorted = fdr_de[order_p]
for i in range(n_g - 2, -1, -1):
    if fdr_sorted[i] > fdr_sorted[i + 1]:
        fdr_sorted[i] = fdr_sorted[i + 1]
fdr_de[order_p] = fdr_sorted

# 3. Calcular Module Membership (kME = correlación de cada transcrito con el Eigengene de su propio módulo en n=9)
samples_9 = list(df_me.index)
X_9 = df_vst.loc[samples_9].values
Z_9 = (X_9 - X_9.mean(axis=0, keepdims=True)) / X_9.std(axis=0, ddof=1, keepdims=True)
ME_std_9 = (df_me.values - df_me.values.mean(axis=0, keepdims=True)) / df_me.values.std(axis=0, ddof=1, keepdims=True)
me_col_to_idx = {col.replace("ME", ""): idx for idx, col in enumerate(df_me.columns)}

kME_own = np.zeros(n_g, dtype=np.float64)
mod_colors = df_mod["Module_Sin_B19"].values
for m_color, m_idx in me_col_to_idx.items():
    mask = mod_colors == m_color
    if np.any(mask):
        kME_own[mask] = (Z_9[:, mask].T @ ME_std_9[:, m_idx]) / float(len(samples_9) - 1)

# 4. Calcular Gene Significance (GS = correlación del transcrito con alcaloides clave en n=8 muestras homologadas)
samples_8 = list(df_meta.index)
X_8 = df_vst.loc[samples_8].values
Z_8 = (X_8 - X_8.mean(axis=0, keepdims=True)) / np.where(X_8.std(axis=0, ddof=1, keepdims=True) == 0, 1.0, X_8.std(axis=0, ddof=1, keepdims=True))

def calc_gs(comp_name):
    y = np.log10(df_meta[comp_name].values)
    y_std = (y - y.mean()) / y.std(ddof=1)
    return (Z_8.T @ y_std) / float(len(samples_8) - 1)

GS_X70T = calc_gs("X70T")
GS_X52T = calc_gs("X52T")
GS_X82P = calc_gs("X82P")
GS_X71T = calc_gs("X71T")

# 5. Consolidar tabla maestra de transcritos
master_df = pd.DataFrame({
    "Transcript_ID": df_vst.columns,
    "Module_Sin_B19": mod_colors,
    "kME_Propio_Modulo": np.round(kME_own, 4),
    "Mean_VST_Alta_Montana": np.round(mean_alta, 3),
    "Mean_VST_Tierras_Bajas": np.round(mean_baja, 3),
    "Mean_VST_Datura": np.round(mean_datura, 3),
    "log2FC_Alta_vs_Baja": np.round(log2FC_Alta_vs_Baja, 3),
    "pvalue_Alta_vs_Baja": np.round(p_val_de, 6),
    "FDR_Alta_vs_Baja": np.round(fdr_de, 5),
    "GS_X70T_log10": np.round(GS_X70T, 4),
    "GS_X52T_log10": np.round(GS_X52T, 4),
    "GS_X82P_log10": np.round(GS_X82P, 4),
    "GS_X71T_log10": np.round(GS_X71T, 4),
})

# Resumen de DEGs por módulo (|log2FC| >= 1.5 y p <= 0.01 entre Alta Montaña y Tierras Bajas)
deg_mask = (np.abs(master_df["log2FC_Alta_vs_Baja"]) >= 1.5) & (master_df["pvalue_Alta_vs_Baja"] <= 0.01)
up_alta = (master_df["log2FC_Alta_vs_Baja"] >= 1.5) & (master_df["pvalue_Alta_vs_Baja"] <= 0.01)
up_baja = (master_df["log2FC_Alta_vs_Baja"] <= -1.5) & (master_df["pvalue_Alta_vs_Baja"] <= 0.01)

print("=== 1. RESUMEN DE EXPRESIÓN DIFERENCIAL: ALTA MONTAÑA (Subclado A) vs. TIERRAS BAJAS (Subclado B) ===")
print(f"Total transcritos diferenciales (|log2FC| >= 1.5, p <= 0.01): {int(deg_mask.sum())}")
print(f"   - Sobre-expresados en Alta Montaña (Subclado A, log2FC >= +1.5): {int(up_alta.sum())}")
print(f"   - Sobre-expresados en Tierras Bajas (Subclado B, log2FC <= -1.5): {int(up_baja.sum())}")

print("\nDistribución de genes diferenciales (Alta Montaña vs. Tierras Bajas) por Módulo WGCNA:")
deg_by_mod = master_df[deg_mask].groupby("Module_Sin_B19").agg(
    Total_DEGs=("Transcript_ID", "count"),
    Up_Alta_Montana=("log2FC_Alta_vs_Baja", lambda s: int((s > 0).sum())),
    Up_Tierras_Bajas=("log2FC_Alta_vs_Baja", lambda s: int((s < 0).sum())),
    Mean_log2FC=("log2FC_Alta_vs_Baja", lambda s: round(float(s.mean()), 2))
).sort_values(by="Total_DEGs", ascending=False)
print(deg_by_mod.head(12).to_string())

# Exportar los Top Hub Genes Candidatos (|kME| >= 0.85, |log2FC| >= 1.5, p <= 0.01 y |GS| >= 0.75)
hub_candidates = master_df[
    (master_df["kME_Propio_Modulo"] >= 0.85) &
    (np.abs(master_df["log2FC_Alta_vs_Baja"]) >= 1.5) &
    (master_df["pvalue_Alta_vs_Baja"] <= 0.01) &
    ((np.abs(master_df["GS_X70T_log10"]) >= 0.75) |
     (np.abs(master_df["GS_X52T_log10"]) >= 0.75) |
     (np.abs(master_df["GS_X82P_log10"]) >= 0.75) |
     (np.abs(master_df["GS_X71T_log10"]) >= 0.75))
].sort_values(by="pvalue_Alta_vs_Baja", ascending=True)

hub_candidates.to_csv("WGCNA_SinB19_04_Top_HubGenes_Candidatos_DEGs_TA.csv", index=False)
print(f"\n=== 2. TOP HUB GENES CANDIDATOS (kME >= 0.85 + DEG Alta/Baja + |GS| >= 0.75): {len(hub_candidates)} transcritos ===")
print(hub_candidates.head(20).to_string(index=False))
