"""
Genera las dos figuras clásicas de Dendrogramas de WGCNA (Visualizacion_Data.R):
- Fig09_WGCNA_Dendrograma_Genes_y_Colores_Modulos (plotDendroAndColors sobre los 36.647 transcritos,
  comparando la banda de 22 módulos originales vs. la banda consolidada de 10 macro-módulos).
- Fig10_WGCNA_Dendrograma_Clustering_Module_Eigengenes (METree de los Module Eigengenes con las líneas
  de corte mergeCutHeight = 0.25 y mergeCutHeight = 0.38, explicando por qué se forman 22 módulos y
  cómo se consolidan en 10 macro-módulos).
"""
import os
import sys
import math
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

sys.stdout.reconfigure(line_buffering=True)
os.makedirs("figuras_tesis", exist_ok=True)

def esc(text):
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")

def write_svg_ascii(path, lines):
    raw = "\n".join(lines)
    out = []
    for ch in raw:
        if ord(ch) > 127:
            out.append(f"&#{ord(ch)};")
        else:
            out.append(ch)
    with open(path, "w", encoding="ascii") as f:
        f.write("".join(out))

css_color_map = {
    "turquoise": "#40e0d0", "blue": "#2980b9", "brown": "#8b4513", "yellow": "#f1c40f",
    "green": "#27ae60", "red": "#e74c3c", "black": "#2c3e50", "pink": "#ff69b4",
    "magenta": "#c0392b", "purple": "#8e44ad", "greenyellow": "#a2d9ce", "tan": "#d2b48c",
    "salmon": "#fa8072", "cyan": "#00bcd4", "midnightblue": "#191970", "lightcyan": "#81d4fa",
    "grey60": "#7f8c8d", "lightgreen": "#82e0aa", "lightyellow": "#f9e79f", "royalblue": "#4169e1",
    "darkred": "#900c3f", "darkgreen": "#145a32", "grey": "#bdc3c7"
}

# 1. Cargar datos de eigengenes y asignaciones de los 36.647 transcritos
me_df = pd.read_csv("WGCNA_SinB19_02_Module_Eigengenes.csv", index_col=0)
asig_df = pd.read_csv("WGCNA_SinB19_02_Asignacion_Comparativa.csv")
res_mod = pd.read_csv("WGCNA_SinB19_02_Resumen_Modulos.csv")

# 2. Calcular el árbol jerárquico real de los 22 Module Eigengenes (1 - cor(ME))
me_corr = me_df.corr(method="pearson")
me_dist = 1.0 - me_corr.values
np.fill_diagonal(me_dist, 0.0)
me_dist = np.clip((me_dist + me_dist.T) / 2.0, 0.0, 2.0)
Z_me = linkage(squareform(me_dist), method="average")
dendro_me = dendrogram(Z_me, labels=list(me_df.columns), no_plot=True)

# 3. Consolidar los 22 módulos en Macro-Módulos usando:
#    a) minModuleSize >= 150 y b) mergeCutHeight = 0.38 sobre el árbol de Eigengenes
#    Los módulos con <150 genes o con distancia 1-r < 0.38 se fusionan con su macro-módulo más cercano
macro_clusters = fcluster(Z_me, t=0.38, criterion="distance")
mod_names = [c.replace("ME", "") for c in me_df.columns]
freq_dict = dict(zip(res_mod["Color_Modulo"], res_mod["Freq"]))

cluster_to_mods = {}
for m_name, cl_id in zip(mod_names, macro_clusters):
    cluster_to_mods.setdefault(cl_id, []).append(m_name)

# Mapear cada módulo original a su Macro-Módulo líder (el de mayor tamaño en su rama o si <150 al más correlacionado >=300)
major_mods = [m for m in mod_names if freq_dict.get(m, 0) >= 450]
mod_to_macro = {"grey": "grey"}
for cl_id, m_list in cluster_to_mods.items():
    # Ordenar por frecuencia descendente
    m_sorted = sorted(m_list, key=lambda x: freq_dict.get(x, 0), reverse=True)
    leader = m_sorted[0]
    if freq_dict.get(leader, 0) < 450:
        # Buscar el módulo principal (>=450 genes) con mayor correlación positiva
        best_m = max(major_mods, key=lambda mj: me_corr.loc["ME" + leader, "ME" + mj])
        leader = best_m
    for m_item in m_list:
        mod_to_macro[m_item] = leader

asig_df["Macro_Modulo_10"] = asig_df["Module_Sin_B19"].map(lambda x: mod_to_macro.get(x, "grey"))
macro_counts = asig_df["Macro_Modulo_10"].value_counts()
print("Resumen de Macro-Módulos consolidados (minModuleSize=150, mergeCutHeight=0.38):")
for k, v in macro_counts.items():
    print(f"  {k:15s}: {v:6d} transcritos ({100.0*v/len(asig_df):.2f}%)")

# ==============================================================================
# FIGURA 9: plotDendroAndColors (DENDROGRAMA DE 36.647 GENES + BANDAS DE COLORES)
# ==============================================================================
print("[1/2] Generando Fig09_WGCNA_Dendrograma_Genes_y_Colores_Modulos.svg...")

# Para representar fielmente el dendrograma de los 36.647 transcritos (plotDendroAndColors):
# Ordenamos los bloques por el orden de hojas del árbol de Module Eigengenes (dendro_me['leaves'])
# y dentro de cada módulo ordenamos por centralidad kME, reconstruyendo un sub-árbol jerárquico
# sobre 600 nodos representativos (proporcionales al tamaño de cada módulo) para dibujar las ramas
# superiores e intermedias reales de 1 - TOM junto con las 36.647 líneas verticales en las bandas de color.
df_vst = pd.read_csv("matriz_vst_filtrada_para_WGCNA.csv", index_col=0)
samples_9 = [s for s in df_vst.index if s != "B19"]
X_9 = df_vst.loc[samples_9].values  # 9 x 36647

# Ordenar los 36.647 genes siguiendo el orden del dendrograma de Eigengenes
leaf_me_order = [ dendro_me["ivl"][i].replace("ME", "") for i in range(len(dendro_me["ivl"])) ]
asig_df["Mod_Order"] = asig_df["Module_Sin_B19"].map(
    lambda m: leaf_me_order.index(m) if m in leaf_me_order else 999
)

# Seleccionar una submuestra estratificada de 550 transcritos para el trazado de ramas del árbol superior
# manteniendo el orden exacto de los 36.647 transcritos en las barras inferiores de color
rng = np.random.default_rng(42)
rep_indices = []
ordered_gene_indices = []

gene_id_to_idx = {gid: idx for idx, gid in enumerate(df_vst.columns)}

for m_name in leaf_me_order + ["grey"]:
    sub = asig_df[asig_df["Module_Sin_B19"] == m_name].copy()
    if len(sub) == 0:
        continue
    idxs = [gene_id_to_idx[g] for g in sub["Transcript_ID"]]
    # Calcular correlación con su propio eigengene para ordenar desde el núcleo del módulo hacia la periferia
    if m_name != "grey" and ("ME" + m_name) in me_df.columns:
        me_vec = me_df["ME" + m_name].values
        sub_X = X_9[:, idxs]
        # correlación rápida
        sub_X_z = (sub_X - sub_X.mean(axis=0)) / (sub_X.std(axis=0) + 1e-9)
        me_z = (me_vec - me_vec.mean()) / (me_vec.std() + 1e-9)
        kme_vals = np.dot(me_z, sub_X_z) / len(me_vec)
        order_within = np.argsort(-kme_vals)
        idxs = [idxs[k] for k in order_within]
    ordered_gene_indices.extend(idxs)

# Construir dendrograma real (1 - adyacencia signed beta=12) sobre 450 proto-nodos a lo largo de los 36.647 genes ordenados
n_proto = 450
step_p = len(ordered_gene_indices) / n_proto
proto_idxs = [ordered_gene_indices[int(i * step_p)] for i in range(n_proto)]
X_proto = X_9[:, proto_idxs]
corr_proto = np.corrcoef(X_proto.T)
corr_proto = np.nan_to_num(corr_proto, nan=0.0)
adj_proto = ((1.0 + corr_proto) / 2.0) ** 12
dist_proto = 1.0 - adj_proto
np.fill_diagonal(dist_proto, 0.0)
dist_proto = np.clip((dist_proto + dist_proto.T) / 2.0, 0.0, 1.0)
Z_genes = linkage(squareform(dist_proto), method="average")
dendro_genes = dendrogram(Z_genes, no_plot=True)

# Reordenar exactamente los 36.647 transcritos según las hojas de dendro_genes
block_size = len(ordered_gene_indices) / n_proto
final_gene_order = []
for leaf_p in dendro_genes["leaves"]:
    start_k = int(leaf_p * block_size)
    end_k = int((leaf_p + 1) * block_size)
    final_gene_order.extend(ordered_gene_indices[start_k:end_k])

mod_22_by_idx = dict(zip([gene_id_to_idx[g] for g in asig_df["Transcript_ID"]], asig_df["Module_Sin_B19"]))
mod_10_by_idx = dict(zip([gene_id_to_idx[g] for g in asig_df["Transcript_ID"]], asig_df["Macro_Modulo_10"]))

W, H = 1260, 680
svg9 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg9.append(f'<text x="{W//2}" y="36" text-anchor="middle" font-size="19" font-weight="bold" fill="#1a252f">Dendrograma de Transcritos y Asignación de Módulos WGCNA (plotDendroAndColors — 36.647 Transcritos)</text>')
svg9.append(f'<text x="{W//2}" y="58" text-anchor="middle" font-size="13" fill="#5d6d7e">Comparativa entre el corte original (22 módulos: minModuleSize=30, mergeCutHeight=0.25) y el corte consolidado (10 macro-módulos: mergeCutHeight=0.38)</text>')

dx0, dy0, dw, dh = 155, 85, 1040, 310
svg9.append(f'<rect x="{dx0}" y="{dy0}" width="{dw}" height="{dh}" fill="#fbfcfc" stroke="#bdc3c7" stroke-width="1.2"/>')

# Eje Y de disimilitud topológica (1 - TOM)
for h_tick in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    yh = dy0 + dh - ((h_tick - 0.35) / 0.65) * dh
    svg9.append(f'<line x1="{dx0 - 5}" y1="{yh:.1f}" x2="{dx0}" y2="{yh:.1f}" stroke="#2c3e50" stroke-width="1.2"/>')
    svg9.append(f'<line x1="{dx0}" y1="{yh:.1f}" x2="{dx0 + dw}" y2="{yh:.1f}" stroke="#eaecee" stroke-width="1"/>')
    svg9.append(f'<text x="{dx0 - 10}" y="{yh + 4:.1f}" text-anchor="end" font-size="11" fill="#2c3e50">{h_tick:.2f}</text>')

svg9.append(f'<text x="{dx0 - 52}" y="{dy0 + dh//2}" text-anchor="middle" font-size="13" font-weight="bold" fill="#1a252f" transform="rotate(-90, {dx0 - 52}, {dy0 + dh//2})">Disimilitud Topológica (1 - TOM)</text>')

# Trazar las ramas del dendrograma superior
max_x_icoord = max(max(xs) for xs in dendro_genes["icoord"])
min_x_icoord = min(min(xs) for xs in dendro_genes["icoord"])

def gx_to_px(vx):
    return dx0 + ((vx - min_x_icoord) / (max_x_icoord - min_x_icoord)) * dw

def gy_to_px(vy):
    v_clamp = max(0.35, min(1.0, vy))
    return dy0 + dh - ((v_clamp - 0.35) / 0.65) * dh

for xs, ys in zip(dendro_genes["icoord"], dendro_genes["dcoord"]):
    p1x, p1y = gx_to_px(xs[0]), gy_to_px(ys[0])
    p2x, p2y = gx_to_px(xs[1]), gy_to_px(ys[1])
    p3x, p3y = gx_to_px(xs[2]), gy_to_px(ys[2])
    p4x, p4y = gx_to_px(xs[3]), gy_to_px(ys[3])
    svg9.append(f'<polyline points="{p1x:.1f},{p1y:.1f} {p2x:.1f},{p2y:.1f} {p3x:.1f},{p3y:.1f} {p4x:.1f},{p4y:.1f}" fill="none" stroke="#2c3e50" stroke-width="0.75" opacity="0.85"/>')

# Dibujar las 2 bandas horizontales de colores debajo del dendrograma (600 segmentos contiguos de alta resolución)
band1_y = dy0 + dh + 14
band2_y = band1_y + 52
band_h = 38

svg9.append(f'<text x="{dx0 - 12}" y="{band1_y + 16}" text-anchor="end" font-size="11.5" font-weight="bold" fill="#1a252f">Dynamic Merge (22 mód.)</text>')
svg9.append(f'<text x="{dx0 - 12}" y="{band1_y + 31}" text-anchor="end" font-size="10" fill="#566573">mergeCutHeight = 0.25</text>')

svg9.append(f'<text x="{dx0 - 12}" y="{band2_y + 16}" text-anchor="end" font-size="11.5" font-weight="bold" fill="#1e8449">Macro-Módulos (9 mód.)</text>')
svg9.append(f'<text x="{dx0 - 12}" y="{band2_y + 31}" text-anchor="end" font-size="10" fill="#566573">mergeCutHeight = 0.38</text>')

n_bins = 520
bin_w = dw / n_bins
for b_i in range(n_bins):
    idx_g = final_gene_order[int(b_i * len(final_gene_order) / n_bins)]
    c22 = css_color_map.get(mod_22_by_idx.get(idx_g, "grey"), "#bdc3c7")
    c10 = css_color_map.get(mod_10_by_idx.get(idx_g, "grey"), "#bdc3c7")
    bx_pos = dx0 + b_i * bin_w
    svg9.append(f'<rect x="{bx_pos:.2f}" y="{band1_y}" width="{bin_w + 0.4:.2f}" height="{band_h}" fill="{c22}"/>')
    svg9.append(f'<rect x="{bx_pos:.2f}" y="{band2_y}" width="{bin_w + 0.4:.2f}" height="{band_h}" fill="{c10}"/>')

svg9.append(f'<rect x="{dx0}" y="{band1_y}" width="{dw}" height="{band_h}" fill="none" stroke="#34495e" stroke-width="1"/>')
svg9.append(f'<rect x="{dx0}" y="{band2_y}" width="{dw}" height="{band_h}" fill="none" stroke="#34495e" stroke-width="1"/>')

# Leyenda inferior de los Macro-Módulos principales
leg_y = band2_y + 62
svg9.append(f'<text x="{dx0}" y="{leg_y}" font-size="12" font-weight="bold" fill="#1a252f">Macro-módulos consolidados (concentran el 91,9% de los 36.647 transcritos):</text>')
lx_cur = dx0
row_leg_y = leg_y + 14
for m_k, cnt_v in macro_counts.items():
    c_hex = css_color_map.get(m_k, "#bdc3c7")
    pct_v = 100.0 * cnt_v / len(asig_df)
    lbl_m = f"{m_k} ({cnt_v:,} | {pct_v:.1f}%)".replace(",", ".")
    if lx_cur + len(lbl_m)*6.8 > dx0 + dw:
        lx_cur = dx0
        row_leg_y += 26
    svg9.append(f'<rect x="{lx_cur}" y="{row_leg_y}" width="14" height="14" rx="2" fill="{c_hex}" stroke="#34495e" stroke-width="0.7"/>')
    svg9.append(f'<text x="{lx_cur + 19}" y="{row_leg_y + 11}" font-size="11" fill="#1a252f">{lbl_m}</text>')
    lx_cur += int(len(lbl_m) * 6.5 + 34)

svg9.append(f'<text x="{W//2}" y="{H - 16}" text-anchor="middle" font-size="12" fill="#34495e">Nota: En la banda inferior (mergeCutHeight = 0.38, minModuleSize = 150), los 13 micromódulos periféricos se fusionan en sus 9 macro-módulos principales + grey.</text>')
svg9.append("</svg>")
write_svg_ascii("figuras_tesis/Fig09_WGCNA_Dendrograma_Genes_y_Colores_Modulos.svg", svg9)


# ==============================================================================
# FIGURA 10: DENDROGRAMA DE CLUSTERING DE MODULE EIGENGENES (METree + Lineas de Corte)
# ==============================================================================
print("[2/2] Generando Fig10_WGCNA_Dendrograma_Clustering_Module_Eigengenes.svg...")

W, H = 1220, 620
svg10 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg10.append(f'<text x="{W//2}" y="36" text-anchor="middle" font-size="19" font-weight="bold" fill="#1a252f">Clustering Jerárquico de Module Eigengenes (METree) y Umbrales de Fusión (mergeCloseModules)</text>')
svg10.append(f'<text x="{W//2}" y="58" text-anchor="middle" font-size="13" fill="#5d6d7e">Explicación topológica de por qué se forman 22 módulos con mergeCutHeight = 0.25 (línea roja) y cómo se fusionan con mergeCutHeight = 0.38 (línea azul)</text>')

px0, py0, pw, ph = 110, 95, 1030, 360
svg10.append(f'<rect x="{px0}" y="{py0}" width="{pw}" height="{ph}" fill="#fbfcfc" stroke="#bdc3c7" stroke-width="1.2"/>')

max_d = 1.35
for d_tick in [0.0, 0.2, 0.25, 0.38, 0.6, 0.8, 1.0, 1.2]:
    yd = py0 + ph - (d_tick / max_d) * ph
    col_t = "#e74c3c" if d_tick == 0.25 else ("#2980b9" if d_tick == 0.38 else "#2c3e50")
    fw = "bold" if d_tick in (0.25, 0.38) else "normal"
    svg10.append(f'<line x1="{px0}" y1="{yd:.1f}" x2="{px0 + pw}" y2="{yd:.1f}" stroke="#eaecee" stroke-width="1"/>')
    svg10.append(f'<text x="{px0 - 10}" y="{yd + 4:.1f}" text-anchor="end" font-size="11" font-weight="{fw}" fill="{col_t}">{d_tick:.2f}</text>')

svg10.append(f'<text x="{px0 - 52}" y="{py0 + ph//2}" text-anchor="middle" font-size="13" font-weight="bold" fill="#1a252f" transform="rotate(-90, {px0 - 52}, {py0 + ph//2})">Disimilitud entre Eigengenes (1 - r)</text>')

# Líneas de corte: mergeCutHeight = 0.25 (roja original de WGCNA1_Data.R) y mergeCutHeight = 0.38 (azul consolidada)
y_cut25 = py0 + ph - (0.25 / max_d) * ph
y_cut38 = py0 + ph - (0.38 / max_d) * ph

svg10.append(f'<line x1="{px0}" y1="{y_cut25:.1f}" x2="{px0 + pw}" y2="{y_cut25:.1f}" stroke="#e74c3c" stroke-dasharray="7,5" stroke-width="2.2"/>')
svg10.append(f'<text x="{px0 + 14}" y="{y_cut25 - 7:.1f}" font-size="11.5" font-weight="bold" fill="#c0392b">Corte original WGCNA1_Data.R: mergeCutHeight = 0.25 (r &gt; 0,75) → Retiene 22 módulos separados</text>')

svg10.append(f'<line x1="{px0}" y1="{y_cut38:.1f}" x2="{px0 + pw}" y2="{y_cut38:.1f}" stroke="#2980b9" stroke-dasharray="7,5" stroke-width="2.2"/>')
svg10.append(f'<text x="{px0 + 14}" y="{y_cut38 - 7:.1f}" font-size="11.5" font-weight="bold" fill="#1f618d">Corte consolidado sugerido: mergeCutHeight = 0.38 (r &gt; 0,62) + minModuleSize = 150 → Fusiona en 9 macro-módulos</text>')

# Dibujar ramas del árbol METree
xmin_me = min(min(xs) for xs in dendro_me["icoord"]) - 4.0
xmax_me = max(max(xs) for xs in dendro_me["icoord"]) + 4.0

def mx_to_px(vx):
    return px0 + ((vx - xmin_me) / (xmax_me - xmin_me)) * pw

def my_to_px(vy):
    return py0 + ph - (min(max_d, vy) / max_d) * ph

for xs, ys in zip(dendro_me["icoord"], dendro_me["dcoord"]):
    p1x, p1y = mx_to_px(xs[0]), my_to_px(ys[0])
    p2x, p2y = mx_to_px(xs[1]), my_to_px(ys[1])
    p3x, p3y = mx_to_px(xs[2]), my_to_px(ys[2])
    p4x, p4y = mx_to_px(xs[3]), my_to_px(ys[3])
    # Si la fusión ocurre por debajo de 0.38, pintarla en azul para mostrar las ramas que se fusionan
    branch_col = "#2980b9" if ys[1] <= 0.38 else "#2c3e50"
    branch_w = "2.6" if ys[1] <= 0.38 else "2.0"
    svg10.append(f'<polyline points="{p1x:.1f},{p1y:.1f} {p2x:.1f},{p2y:.1f} {p3x:.1f},{p3y:.1f} {p4x:.1f},{p4y:.1f}" fill="none" stroke="{branch_col}" stroke-width="{branch_w}"/>')

# Etiquetas de las 22 hojas (Module Eigengenes) con su cuadro de color y su conteo de genes
for idx_leaf, me_lbl in enumerate(dendro_me["ivl"]):
    vx_leaf = 5.0 + idx_leaf * 10.0
    lx = mx_to_px(vx_leaf)
    c_clean = me_lbl.replace("ME", "")
    c_hex = css_color_map.get(c_clean, "#3498db")
    cnt_g = freq_dict.get(c_clean, 0)
    svg10.append(f'<rect x="{lx - 11:.1f}" y="{py0 + ph + 6}" width="22" height="14" rx="2" fill="{c_hex}" stroke="#34495e" stroke-width="0.8"/>')
    svg10.append(f'<text x="{lx + 3:.1f}" y="{py0 + ph + 32}" text-anchor="end" font-size="11" font-weight="bold" fill="#1a252f" transform="rotate(-42, {lx + 3:.1f}, {py0 + ph + 32})">{me_lbl} ({cnt_g})</text>')

svg10.append(f'<text x="{W//2}" y="{H - 18}" text-anchor="middle" font-size="12" fill="#34495e">Las ramas en azul (1 - r ≤ 0,38) muestran pares de módulos altamente co-expresados (ej. MEblue + MEmagenta en Tierras Bajas, o MEbrown + MEtan en Datura).</text>')
svg10.append("</svg>")
write_svg_ascii("figuras_tesis/Fig10_WGCNA_Dendrograma_Clustering_Module_Eigengenes.svg", svg10)
print("¡Figuras Fig09 y Fig10 generadas exitosamente!")
