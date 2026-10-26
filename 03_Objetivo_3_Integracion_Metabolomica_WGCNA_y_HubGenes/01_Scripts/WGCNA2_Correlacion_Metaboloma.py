"""
Fase 2 WGCNA: Cruce de Module Eigengenes (22 módulos sin B19) contra el Metaboloma Real (Metaboloma_WGNA_TESIS.xlsx)
Tesis Hologenómica Clado Datura — Módulo 4
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)

# 1. Cargar Module Eigengenes (n=9 sin B19: B9, B10, B11, B12, B14, B15, B17, B21, B27)
me_df = pd.read_csv("WGCNA_SinB19_02_Module_Eigengenes.csv", index_col=0)
# Cargar el metaboloma corregido fiel a Metaboloma_WGNA_TESIS.xlsx
meta_df = pd.read_csv("fenotipos_TA_corregidos_desde_Excel.csv", index_col=0)

print("=== 1. HOMOLOGACIÓN DE MUESTRAS EXPRESIÓN vs. METABOLÓMICA ===")
print("Muestras en Expresión (sin B19):", list(me_df.index))
print("Muestras en Metabolómica       :", list(meta_df.index))
solo_expr = set(me_df.index) - set(meta_df.index)
solo_meta = set(meta_df.index) - set(me_df.index)
print(f"Exclusivas de Expresión   : {sorted(solo_expr)} (Datura wrightii B27)")
print(f"Exclusivas de Metabolómica: {sorted(solo_meta)} (Outgroups B19 Atropa y B20 Ipomoea)")

common_samples = [s for s in me_df.index if s in meta_df.index]
print(f"Muestras comunes homologadas (n={len(common_samples)}): {common_samples}")

me_aligned = me_df.loc[common_samples].copy()
meta_aligned = meta_df.loc[common_samples].copy()
assert list(me_aligned.index) == list(meta_aligned.index), "Error: las filas no coinciden"

# Exportar matriz metabolómica limpia, corregida y homologada para R
meta_df.to_csv("fenotipos_TA_corregidos_para_WGCNA.csv", index=True)
meta_aligned.to_csv("fenotipos_TA_homologados_n8_SinOutgroups.csv", index=True)

# Separar los 76 Alcaloides Tropánicos (X6T..X81T) y los 3 Precursores (X82P..X84P)
ta_cols = [c for c in meta_aligned.columns if c.endswith("T")]
prec_cols = [c for c in meta_aligned.columns if c.endswith("P")]
print(f"\nCompuestos en la matriz: {len(ta_cols)} Alcaloides Tropánicos (T) + {len(prec_cols)} Precursores (P) = {meta_aligned.shape[1]} total")

# 2. Detectar clusters de picos altamente correlacionados (r > 0.85, posibles isómeros / fragmentos en fuente)
corr_meta = meta_aligned[ta_cols].corr()
high_corr_pairs = []
for i in range(len(ta_cols)):
    for j in range(i + 1, len(ta_cols)):
        r_ij = corr_meta.iloc[i, j]
        if r_ij > 0.85:
            high_corr_pairs.append((ta_cols[i], ta_cols[j], round(float(r_ij), 4)))

high_corr_pairs.sort(key=lambda x: -x[2])
print(f"\n=== 2. PARES DE ALCALOIDES CON r > 0.85 (Candidatos a Fragmentos en Fuente / Isómeros): {len(high_corr_pairs)} pares ===")
for p1, p2, r_v in high_corr_pairs[:15]:
    print(f"   {p1:5s} <-> {p2:5s} : r = {r_v:.4f}")

# 3. Calcular correlación de Pearson (r) y p-valor (t-Student) entre los 22 Módulos y los 79 Compuestos
# Evaluamos sobre log10(abundancia) (estándar metabolómica MS para estabilizar asimetría de áreas de pico) y escala lineal
meta_log10 = np.log10(meta_aligned)

n_s = len(common_samples)
r_mat = pd.DataFrame(index=me_aligned.columns, columns=meta_aligned.columns, dtype=float)
p_mat = pd.DataFrame(index=me_aligned.columns, columns=meta_aligned.columns, dtype=float)

top_hits = []
for mod in me_aligned.columns:
    x_vec = me_aligned[mod].values
    for comp in meta_aligned.columns:
        y_vec = meta_log10[comp].values
        r_val, p_val = stats.pearsonr(x_vec, y_vec)
        r_mat.loc[mod, comp] = round(float(r_val), 4)
        p_mat.loc[mod, comp] = round(float(p_val), 5)
        if abs(r_val) >= 0.70 and p_val <= 0.05:
            top_hits.append({
                "Modulo": mod,
                "Compuesto": comp,
                "Tipo": "Alcaloide_TA" if comp.endswith("T") else "Precursor",
                "r_Pearson_log10": round(float(r_val), 4),
                "p_valor": round(float(p_val), 5),
                "r_Pearson_lineal": round(float(stats.pearsonr(x_vec, meta_aligned[comp].values)[0]), 4)
            })

r_mat.to_csv("WGCNA_SinB19_03_Correlacion_Modulos_Metabolitos.csv")
p_mat.to_csv("WGCNA_SinB19_03_Pvalores_Modulos_Metabolitos.csv")

hits_df = pd.DataFrame(top_hits).sort_values(by="p_valor", ascending=True)
hits_df.to_csv("WGCNA_SinB19_03_Top_Asociaciones_Significativas.csv", index=False)

print(f"\n=== 3. ASOCIACIONES SIGNIFICATIVAS MÓDULO-ALCALOIDE (|r| >= 0.70, p <= 0.05): {len(hits_df)} asociaciones ===")
print(hits_df.head(25).to_string(index=False))

print("\nResumen de asociaciones significativas (|r| >= 0.70, p <= 0.05) por Módulo:")
print(hits_df["Modulo"].value_counts().to_string())
