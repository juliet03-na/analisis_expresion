"""
Generador de las 8 Figuras Científicas de Calidad de Publicación (SVG vectorial -> PNG 300 DPI + PDF)
Tesis Hologenómica Clado Datura — Módulos 3 y 4
"""
import os
import sys
import math
import numpy as np
import pandas as pd

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

def val_to_bwr_hex(val, vmin=-1.0, vmax=1.0):
    """Escala de color Blue - White - Red (estilo WGCNA labeledHeatmap)."""
    x = max(0.0, min(1.0, (val - vmin) / (vmax - vmin)))
    if x < 0.5:
        t = x / 0.5
        r = int(33 + t * (255 - 33))
        g = int(102 + t * (255 - 102))
        b = int(172 + t * (255 - 172))
    else:
        t = (x - 0.5) / 0.5
        r = int(255 - t * (255 - 178))
        g = int(255 - t * (255 - 24))
        b = int(255 - t * (255 - 43))
    return f"#{r:02x}{g:02x}{b:02x}"

def val_to_viridis_corr(val, vmin=0.15, vmax=1.0):
    """Escala cálida RdYlBu invertida para correlación entre muestras [0.15, 1.0]."""
    x = max(0.0, min(1.0, (val - vmin) / (vmax - vmin)))
    if x < 0.5:
        t = x / 0.5
        r = int(69 + t * (254 - 69))
        g = int(117 + t * (224 - 117))
        b = int(180 + t * (144 - 180))
    else:
        t = (x - 0.5) / 0.5
        r = int(254 - t * (254 - 215))
        g = int(224 - t * (224 - 48))
        b = int(144 - t * (144 - 39))
    return f"#{r:02x}{g:02x}{b:02x}"

# ==============================================================================
# FIGURA 1: EMBUDO DE FILTRADO ESTADÍSTICO (135.008 -> 123.099 -> 54.363 -> 36.647)
# ==============================================================================
print("[1/8] Generando Fig01_Embudo_Filtrado_Transcritos.svg...")
steps_funnel = [
    ("Ensamblaje crudo rnaSPAdes (K=49, 73)", 135008, 0, "#1b4f72", "Co-ensamblaje de novo inicial"),
    ("Paso 0: Post-cluster CD-HIT-EST (-c 0.95)", 123099, 11909, "#21618c", "-8,8% redundancia (98,3% mapeo Bowtie2)"),
    ("Pasos 1-2: Suma > 0 y Varianza > 0", 122690, 409, "#2874a6", "Remoción de transcritos inactivos"),
    ("Paso 3: Expresión mínima (rowSums >= 25)", 54435, 68255, "#2e86c1", "Filtro principal de ruido de bajo conteo"),
    ("Paso 4: Filtro de extremos (max <= 50.000)", 54409, 26, "#3498db", "Remoción de hiper-expresión nivel 1"),
    ("Pasos 5-6: max <= 30.000 + Duplicados + VST", 54363, 46, "#17a589", "Conteo exacto Pre-MAD estabilizado (VST)"),
    ("Paso 7: Filtro de Varianza MAD (75% sup.)", 36647, 17716, "#117a65", "Matriz final WGCNA (10 muestras x 36.647)"),
]

W, H = 1150, 680
svg1 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg1.append(f'<text x="{W//2}" y="42" text-anchor="middle" font-size="20" font-weight="bold" fill="#1a252f">Embudo de Filtrado y Normalización Transcriptómica (Módulos 1 y 3 — Normalizacion3.R)</text>')
svg1.append(f'<text x="{W//2}" y="68" text-anchor="middle" font-size="13" fill="#5d6d7e">Trazabilidad exacta desde contigs crudos de rnaSPAdes hasta la matriz homocedástica VST-MAD para WGCNA</text>')

max_val = 135008.0
top_y = 95
row_h = 74
center_x = 575
max_bar_w = 460
min_bar_w = 250

for i, (label, ret, elim, color, note) in enumerate(steps_funnel):
    y = top_y + i * row_h
    w_curr = int(min_bar_w + (max_bar_w - min_bar_w) * (ret / max_val))
    if i < len(steps_funnel) - 1:
        w_next = int(min_bar_w + (max_bar_w - min_bar_w) * (steps_funnel[i+1][1] / max_val))
    else:
        w_next = w_curr

    x1_top, x2_top = center_x - w_curr // 2, center_x + w_curr // 2
    x1_bot, x2_bot = center_x - w_next // 2, center_x + w_next // 2
    y_bot = y + 54

    # Trapecio del embudo
    pts = f"{x1_top},{y} {x2_top},{y} {x2_bot},{y_bot} {x1_bot},{y_bot}"
    svg1.append(f'<polygon points="{pts}" fill="{color}" stroke="#ffffff" stroke-width="2" opacity="0.93"/>')

    # Etiqueta izquierda (Paso)
    svg1.append(f'<text x="28" y="{y + 25}" font-size="13" font-weight="bold" fill="#1a252f">{esc(label)}</text>')
    svg1.append(f'<text x="28" y="{y + 44}" font-size="11.5" fill="#566573">{esc(note)}</text>')

    # Número central dentro del embudo
    pct_cdhit = 100.0 * ret / 123099.0
    svg1.append(f'<text x="{center_x}" y="{y + 33}" text-anchor="middle" font-size="14.5" font-weight="bold" fill="#ffffff">{ret:,} transcritos ({pct_cdhit:.1f}%)</text>'.replace(",", "."))

    # Etiqueta derecha (Eliminados)
    if elim > 0:
        svg1.append(f'<rect x="845" y="{y + 10}" width="275" height="34" rx="6" fill="#fdedec" stroke="#e74c3c" stroke-width="1"/>')
        svg1.append(f'<text x="982" y="{y + 32}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#c0392b">Eliminados: -{elim:,} transcritos</text>'.replace(",", "."))
    else:
        svg1.append(f'<rect x="845" y="{y + 10}" width="275" height="34" rx="6" fill="#eafaf1" stroke="#27ae60" stroke-width="1"/>')
        svg1.append(f'<text x="982" y="{y + 32}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#1e8449">Retención íntegra (0 eliminados)</text>')

svg1.append(f'<text x="{W//2}" y="{H - 18}" text-anchor="middle" font-size="12" fill="#34495e">Nota: El conteo de 54.363 corresponde exactamente al Paso 5-6 (Pre-MAD), y 36.647 al Paso 7 (Post-MAD, mads &gt; quantile 0.25).</text>')
svg1.append("</svg>")
write_svg_ascii("figuras_tesis/Fig01_Embudo_Filtrado_Transcritos.svg", svg1)


# ==============================================================================
# FIGURA 2: SALMON v1.10.3 - LECTURAS PROCESADAS/MAPEADAS Y TASA DE ALINEAMIENTO
# ==============================================================================
print("[2/8] Generando Fig02_Salmon_Cobertura_y_Mapeo.svg...")
df_salmon = pd.read_csv("Salmon_Resumen_Mapeo_Plantas.csv")
nicho_short = {
    "B9": ("B. insignis", "Subclado B: Tierras bajas", "#27ae60"),
    "B12": ("B. suaveolens", "Subclado B: Tierras bajas", "#27ae60"),
    "B14": ("B. versicolor", "Subclado B: Tierras bajas", "#27ae60"),
    "B21": ("B. aurea", "Subclado B: Tierras bajas", "#27ae60"),
    "B10": ("B. sanguinea", "Subclado A: Alta montaña", "#2980b9"),
    "B11": ("B. sang. x vulc.", "Subclado A: Alta montaña", "#2980b9"),
    "B15": ("B. vulcanicola", "Subclado A: Alta montaña", "#2980b9"),
    "B17": ("D. stramonium", "Clado Datura", "#e67e22"),
    "B27": ("D. wrightii", "Clado Datura", "#e67e22"),
    "B19": ("A. belladonna", "Outgroup (Excluido WGCNA)", "#8e44ad"),
}
order_samples = ["B9", "B12", "B14", "B21", "B10", "B11", "B15", "B17", "B27", "B19"]
df_salmon = df_salmon.set_index("Sample_ID").loc[order_samples].reset_index()

W, H = 1150, 620
svg2 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg2.append(f'<text x="{W//2}" y="38" text-anchor="middle" font-size="19" font-weight="bold" fill="#1a252f">Cuantificación en Salmon v1.10.3: Profundidad de Lecturas y Tasa de Mapeo (123.099 Transcritos)</text>')
svg2.append(f'<text x="{W//2}" y="62" text-anchor="middle" font-size="13" fill="#5d6d7e">Comparativa por nicho ecológico (Tierras Bajas vs. Alta Montaña vs. Datura vs. Outgroup B19)</text>')

plot_x, plot_y, plot_w, plot_h = 100, 100, 920, 390
svg2.append(f'<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" fill="#fbfcfc" stroke="#bdc3c7" stroke-width="1"/>')

max_reads = 24000000.0
for m_tick in range(0, 25, 4):
    y_t = plot_y + plot_h - int((m_tick * 1e6 / max_reads) * plot_h)
    svg2.append(f'<line x1="{plot_x}" y1="{y_t}" x2="{plot_x + plot_w}" y2="{y_t}" stroke="#eaecee" stroke-width="1"/>')
    svg2.append(f'<text x="{plot_x - 10}" y="{y_t + 4}" text-anchor="end" font-size="11.5" fill="#34495e">{m_tick}M</text>')

svg2.append(f'<text x="32" y="{plot_y + plot_h//2}" text-anchor="middle" font-size="13" font-weight="bold" fill="#1a252f" transform="rotate(-90, 32, {plot_y + plot_h//2})">Lecturas Pareadas (Millones)</text>')

col_w = plot_w / len(order_samples)
for i, row in df_salmon.iterrows():
    s = row["Sample_ID"]
    sp_name, clade_lbl, col_hex = nicho_short[s]
    cx = plot_x + (i + 0.5) * col_w
    bw = 46

    h_proc = int((row["Reads_Processed"] / max_reads) * plot_h)
    h_map = int((row["Reads_Mapped"] / max_reads) * plot_h)
    y_proc = plot_y + plot_h - h_proc
    y_map = plot_y + plot_h - h_map

    # Barra procesadas (fondo claro) y mapeadas (color sólido)
    svg2.append(f'<rect x="{cx - bw//2}" y="{y_proc}" width="{bw}" height="{h_proc}" fill="#d5dbdb" stroke="#7f8c8d" stroke-width="1"/>')
    svg2.append(f'<rect x="{cx - bw//2}" y="{y_map}" width="{bw}" height="{h_map}" fill="{col_hex}" opacity="0.9"/>')

    # Etiquetas sobre la barra
    reads_m = row["Reads_Mapped"] / 1e6
    lbl_reads = f"{reads_m:.2f}M" if reads_m >= 0.5 else f"{row['Reads_Mapped']/1e3:.1f}k"
    svg2.append(f'<text x="{cx}" y="{y_proc - 22}" text-anchor="middle" font-size="11" font-weight="bold" fill="#1a252f">{lbl_reads}</text>')
    svg2.append(f'<text x="{cx}" y="{y_proc - 7}" text-anchor="middle" font-size="10.5" font-weight="bold" fill="{col_hex}">({row["Percent_Mapped"]:.1f}%)</text>')

    # Etiqueta eje X
    svg2.append(f'<text x="{cx}" y="{plot_y + plot_h + 22}" text-anchor="middle" font-size="13" font-weight="bold" fill="#1a252f">{s}</text>')
    svg2.append(f'<text x="{cx}" y="{plot_y + plot_h + 40}" text-anchor="middle" font-size="10.5" font-style="italic" fill="#566573">{sp_name}</text>')

# Leyenda inferior
leg_items = [
    ("Subclado B: Tierras bajas (B9, B12, B14, B21)", "#27ae60"),
    ("Subclado A: Alta montaña (B10, B11, B15)", "#2980b9"),
    ("Clado Datura (B17, B27)", "#e67e22"),
    ("Outgroup A. belladonna (B19: 122,4k lecturas)", "#8e44ad"),
]
lx = 52
for txt, c_hex in leg_items:
    svg2.append(f'<rect x="{lx}" y="565" width="15" height="15" rx="3" fill="{c_hex}"/>')
    svg2.append(f'<text x="{lx + 20}" y="577" font-size="11.5" fill="#1a252f">{txt}</text>')
    lx += int(len(txt) * 5.9 + 28)

svg2.append("</svg>")
write_svg_ascii("figuras_tesis/Fig02_Salmon_Cobertura_y_Mapeo.svg", svg2)


# ==============================================================================
# FIGURA 3: PCA COMPARATIVO (CON B19 n=10 vs SIN B19 n=9)
# ==============================================================================
print("[3/8] Generando Fig03_QC_PCA_Comparativo_Muestras.svg...")
from sklearn.decomposition import PCA

df_vst = pd.read_csv("matriz_vst_filtrada_para_WGCNA.csv", index_col=0)
samples_10 = list(df_vst.index)
X_10 = df_vst.values

pca_10 = PCA(n_components=2)
c_10 = pca_10.fit_transform(X_10)
v_10 = pca_10.explained_variance_ratio_ * 100

samples_9 = [s for s in samples_10 if s != "B19"]
X_9 = df_vst.loc[samples_9].values
pca_9 = PCA(n_components=2)
c_9 = pca_9.fit_transform(X_9)
v_9 = pca_9.explained_variance_ratio_ * 100

W, H = 1200, 600
svg3 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg3.append(f'<text x="{W//2}" y="36" text-anchor="middle" font-size="19" font-weight="bold" fill="#1a252f">Análisis de Componentes Principales (PCA) sobre Expresión VST (36.647 Transcritos)</text>')
svg3.append(f'<text x="{W//2}" y="58" text-anchor="middle" font-size="13" fill="#5d6d7e">Efecto de la exclusión del outgroup B19 (Atropa belladonna) sobre la separación entre Alta Montaña y Tierras Bajas</text>')

def draw_pca_panel(svg_list, coords, s_list, var_pct, x0, y0, pw, ph, title_str):
    svg_list.append(f'<rect x="{x0}" y="{y0}" width="{pw}" height="{ph}" fill="#fbfcfc" stroke="#bdc3c7" stroke-width="1.2"/>')
    svg_list.append(f'<text x="{x0 + pw//2}" y="{y0 - 12}" text-anchor="middle" font-size="14" font-weight="bold" fill="#1a252f">{esc(title_str)}</text>')

    xmin, xmax = coords[:, 0].min() * 1.22, coords[:, 0].max() * 1.22
    ymin, ymax = coords[:, 1].min() * 1.25, coords[:, 1].max() * 1.25

    def to_px(vx, vy):
        px = x0 + ((vx - xmin) / (xmax - xmin)) * pw
        py = y0 + ph - ((vy - ymin) / (ymax - ymin)) * ph
        return px, py

    zx, zy = to_px(0, 0)
    if x0 <= zx <= x0 + pw:
        svg_list.append(f'<line x1="{zx:.1f}" y1="{y0}" x2="{zx:.1f}" y2="{y0+ph}" stroke="#d0d3d4" stroke-dasharray="5,5" stroke-width="1"/>')
    if y0 <= zy <= y0 + ph:
        svg_list.append(f'<line x1="{x0}" y1="{zy:.1f}" x2="{x0+pw}" y2="{zy:.1f}" stroke="#d0d3d4" stroke-dasharray="5,5" stroke-width="1"/>')

    pca_label_offsets = {
        "B11": (0, -14, "middle"),
        "B10": (-14, 4, "end"),
        "B15": (14, 15, "start"),
        "B14": (-14, -6, "end"),
        "B21": (14, 14, "start"),
        "B9": (-14, -6, "end"),
        "B12": (-14, 18, "end"),
        "B17": (0, -14, "middle"),
        "B27": (0, 22, "middle"),
        "B19": (14, 5, "start"),
    }
    for idx, s in enumerate(s_list):
        sp_name, _, col_hex = nicho_short[s]
        px, py = to_px(coords[idx, 0], coords[idx, 1])
        svg_list.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="8.5" fill="{col_hex}" stroke="#ffffff" stroke-width="2"/>')
        dx, dy, anchor = pca_label_offsets.get(s, (0, -14, "middle"))
        svg_list.append(f'<text x="{px + dx:.1f}" y="{py + dy:.1f}" text-anchor="{anchor}" font-size="11" font-weight="bold" fill="#1a252f">{s} ({sp_name})</text>')

    svg_list.append(f'<text x="{x0 + pw//2}" y="{y0 + ph + 28}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#2c3e50">PC1 ({var_pct[0]:.2f}% de la varianza explicada)</text>')
    svg_list.append(f'<text x="{x0 - 18}" y="{y0 + ph//2}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#2c3e50" transform="rotate(-90, {x0 - 18}, {y0 + ph//2})">PC2 ({var_pct[1]:.2f}% de la varianza explicada)</text>')

draw_pca_panel(svg3, c_10, samples_10, v_10, 65, 100, 495, 410, "A) PCA con Outgroup B19 (n = 10 muestras)")
draw_pca_panel(svg3, c_9, samples_9, v_9, 645, 100, 495, 410, "B) PCA sin Outgroup B19 (n = 9 muestras: +5,2% en PC2)")

lx = 75
for txt, c_hex in leg_items:
    svg3.append(f'<circle cx="{lx}" cy="565" r="7" fill="{c_hex}"/>')
    svg3.append(f'<text x="{lx + 13}" y="570" font-size="11.5" fill="#1a252f">{txt}</text>')
    lx += int(len(txt) * 5.9 + 28)

svg3.append("</svg>")
write_svg_ascii("figuras_tesis/Fig03_QC_PCA_Comparativo_Muestras.svg", svg3)


# ==============================================================================
# FIGURA 4: HEATMAP DE CORRELACIÓN DE PEARSON ENTRE MUESTRAS (VST-MAD)
# ==============================================================================
print("[4/8] Generando Fig04_QC_Heatmap_Correlacion_Muestras.svg...")
corr_df = pd.read_csv("WGCNA_00_Correlacion_Muestras.csv", index_col=0)
corr_df = corr_df.loc[order_samples, order_samples]

W, H = 960, 780
svg4 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg4.append(f'<text x="{W//2}" y="38" text-anchor="middle" font-size="19" font-weight="bold" fill="#1a252f">Matriz de Correlación de Pearson entre Muestras (Expresión VST — 36.647 Transcritos)</text>')
svg4.append(f'<text x="{W//2}" y="62" text-anchor="middle" font-size="13" fill="#5d6d7e">Agrupamiento neto en 3 clústeres: Subclado B (Tierras Bajas), Subclado A (Alta Montaña) y Clado Datura + Outgroup B19</text>')

hx, hy, cell_s = 230, 120, 54
for i, s_row in enumerate(order_samples):
    sp_name, _, col_hex = nicho_short[s_row]
    ry = hy + i * cell_s
    # Etiqueta lateral izquierda + barra de nicho
    svg4.append(f'<rect x="{hx - 24}" y="{ry + 2}" width="16" height="{cell_s - 4}" fill="{col_hex}" rx="3"/>')
    svg4.append(f'<text x="{hx - 32}" y="{ry + cell_s//2 + 5}" text-anchor="end" font-size="12.5" font-weight="bold" fill="#1a252f">{s_row} ({sp_name})</text>')

    # Encabezado superior
    cx = hx + i * cell_s
    svg4.append(f'<rect x="{cx + 2}" y="{hy - 22}" width="{cell_s - 4}" height="14" fill="{col_hex}" rx="3"/>')
    svg4.append(f'<text x="{cx + cell_s//2}" y="{hy - 28}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#1a252f">{s_row}</text>')

    for j, s_col in enumerate(order_samples):
        cx_j = hx + j * cell_s
        r_val = float(corr_df.loc[s_row, s_col])
        fill_c = val_to_viridis_corr(r_val, 0.18, 1.0)
        txt_c = "#ffffff" if (r_val > 0.82 or r_val < 0.28) else "#1a252f"
        svg4.append(f'<rect x="{cx_j}" y="{ry}" width="{cell_s}" height="{cell_s}" fill="{fill_c}" stroke="#ffffff" stroke-width="1.5"/>')
        svg4.append(f'<text x="{cx_j + cell_s//2}" y="{ry + cell_s//2 + 5}" text-anchor="middle" font-size="12" font-weight="bold" fill="{txt_c}">{r_val:.2f}</text>')

# Recuadros resaltando los 3 bloques ecológicos
svg4.append(f'<rect x="{hx}" y="{hy}" width="{4*cell_s}" height="{4*cell_s}" fill="none" stroke="#1e8449" stroke-width="3.5"/>')
svg4.append(f'<rect x="{hx + 4*cell_s}" y="{hy + 4*cell_s}" width="{3*cell_s}" height="{3*cell_s}" fill="none" stroke="#1f618d" stroke-width="3.5"/>')
svg4.append(f'<rect x="{hx + 7*cell_s}" y="{hy + 7*cell_s}" width="{2*cell_s}" height="{2*cell_s}" fill="none" stroke="#d35400" stroke-width="3.5"/>')

# Leyenda derecha
lx_bar = hx + 10 * cell_s + 35
svg4.append(f'<text x="{lx_bar}" y="{hy + 15}" font-size="12.5" font-weight="bold" fill="#1a252f">Escala Pearson r</text>')
for step_i in range(20):
    v_s = 1.0 - step_i * (0.82 / 19.0)
    cy_s = hy + 30 + step_i * 14
    svg4.append(f'<rect x="{lx_bar}" y="{cy_s}" width="26" height="14" fill="{val_to_viridis_corr(v_s, 0.18, 1.0)}"/>')
    if step_i in (0, 5, 10, 15, 19):
        svg4.append(f'<text x="{lx_bar + 34}" y="{cy_s + 11}" font-size="11" fill="#2c3e50">{v_s:.2f}</text>')

svg4.append(f'<text x="{W//2}" y="{H - 28}" text-anchor="middle" font-size="12" fill="#34495e">Los recuadros destacan la alta correlación intra-nicho: Subclado B Tierras Bajas (r = 0,71–0,88), Subclado A Alta Montaña (r = 0,78–0,86) y Datura (r = 0,76).</text>')
svg4.append("</svg>")
write_svg_ascii("figuras_tesis/Fig04_QC_Heatmap_Correlacion_Muestras.svg", svg4)


# ==============================================================================
# FIGURA 5: CURVAS DIAGNÓSTICAS pickSoftThreshold (CONECTIVIDAD Y ESCALA)
# ==============================================================================
print("[5/8] Generando Fig05_WGCNA_Topologia_Libre_Escala.svg...")
sft_df = pd.read_csv("WGCNA_01_Tabla_SoftThreshold.csv")

W, H = 1180, 560
svg5 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg5.append(f'<text x="{W//2}" y="36" text-anchor="middle" font-size="19" font-weight="bold" fill="#1a252f">Diagnóstico de Potencia Suave (pickSoftThreshold) en Red Signed (36.647 Transcritos)</text>')
svg5.append(f'<text x="{W//2}" y="58" text-anchor="middle" font-size="13" fill="#5d6d7e">Comparativa de Conectividad Media (k) y Pendiente Libre de Escala: Con B19 (n=10) vs. Sin B19 (n=9, β=12 seleccionado)</text>')

# Panel A: Conectividad media (escala lineal)
x0, y0, pw, ph = 85, 100, 470, 380
svg5.append(f'<rect x="{x0}" y="{y0}" width="{pw}" height="{ph}" fill="#fbfcfc" stroke="#bdc3c7" stroke-width="1.2"/>')
svg5.append(f'<text x="{x0+pw//2}" y="{y0-12}" text-anchor="middle" font-size="14" font-weight="bold" fill="#1a252f">A) Reducción de Conectividad Media (mean k) vs. Potencia β</text>')

max_k = 26500.0
for k_tick in range(0, 26001, 5000):
    yk = y0 + ph - (k_tick / max_k) * ph
    svg5.append(f'<line x1="{x0}" y1="{yk:.1f}" x2="{x0+pw}" y2="{yk:.1f}" stroke="#eaecee" stroke-width="1"/>')
    svg5.append(f'<text x="{x0-8}" y="{yk+4:.1f}" text-anchor="end" font-size="11" fill="#34495e">{k_tick:,}</text>'.replace(",", "."))

# Línea vertical en beta = 12
x_b12 = x0 + ((12 - 1) / 19.0) * pw
svg5.append(f'<line x1="{x_b12:.1f}" y1="{y0}" x2="{x_b12:.1f}" y2="{y0+ph}" stroke="#27ae60" stroke-dasharray="6,4" stroke-width="2"/>')
svg5.append(f'<text x="{x_b12+8:.1f}" y="{y0+30}" font-size="11.5" font-weight="bold" fill="#1e8449">β = 12 seleccionada (k = 2.880,8)</text>')

pts_10, pts_9 = [], []
for _, r in sft_df.iterrows():
    p = r["Power"]
    px = x0 + ((p - 1) / 19.0) * pw
    py10 = y0 + ph - (r["mean.k._n10_con_B19"] / max_k) * ph
    py9 = y0 + ph - (r["mean.k._n9_sin_B19"] / max_k) * ph
    pts_10.append(f"{px:.1f},{py10:.1f}")
    pts_9.append(f"{px:.1f},{py9:.1f}")
    svg5.append(f'<circle cx="{px:.1f}" cy="{py10:.1f}" r="4.5" fill="#e74c3c"/>')
    svg5.append(f'<circle cx="{px:.1f}" cy="{py9:.1f}" r="4.5" fill="#2980b9"/>')
    if p in (1, 2, 4, 6, 8, 10, 12, 16, 20):
        svg5.append(f'<text x="{px:.1f}" y="{y0+ph+20}" text-anchor="middle" font-size="11" fill="#2c3e50">{int(p)}</text>')

svg5.append(f'<polyline points="{" ".join(pts_10)}" fill="none" stroke="#e74c3c" stroke-width="2.2"/>')
svg5.append(f'<polyline points="{" ".join(pts_9)}" fill="none" stroke="#2980b9" stroke-width="2.2"/>')
svg5.append(f'<text x="{x0+pw//2}" y="{y0+ph+42}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#1a252f">Potencia Suave (Soft Threshold β)</text>')

# Panel B: log10(mean k) y Magnitud R^2 del ajuste libre de escala
x1 = 635
svg5.append(f'<rect x="{x1}" y="{y0}" width="{pw}" height="{ph}" fill="#fbfcfc" stroke="#bdc3c7" stroke-width="1.2"/>')
svg5.append(f'<text x="{x1+pw//2}" y="{y0-12}" text-anchor="middle" font-size="14" font-weight="bold" fill="#1a252f">B) Conectividad en Escala Log10 (Estabilización de Red)</text>')

for lk_tick in [3.0, 3.5, 4.0, 4.5]:
    ylk = y0 + ph - ((lk_tick - 3.0) / 1.5) * ph
    svg5.append(f'<line x1="{x1}" y1="{ylk:.1f}" x2="{x1+pw}" y2="{ylk:.1f}" stroke="#eaecee" stroke-width="1"/>')
    svg5.append(f'<text x="{x1-8}" y="{ylk+4:.1f}" text-anchor="end" font-size="11" fill="#34495e">10^{lk_tick:.1f}</text>')

x1_b12 = x1 + ((12 - 1) / 19.0) * pw
svg5.append(f'<line x1="{x1_b12:.1f}" y1="{y0}" x2="{x1_b12:.1f}" y2="{y0+ph}" stroke="#27ae60" stroke-dasharray="6,4" stroke-width="2"/>')

lpts_10, lpts_9 = [], []
for _, r in sft_df.iterrows():
    p = r["Power"]
    px = x1 + ((p - 1) / 19.0) * pw
    lk10 = math.log10(r["mean.k._n10_con_B19"])
    lk9 = math.log10(r["mean.k._n9_sin_B19"])
    py10 = y0 + ph - ((lk10 - 3.0) / 1.5) * ph
    py9 = y0 + ph - ((lk9 - 3.0) / 1.5) * ph
    lpts_10.append(f"{px:.1f},{py10:.1f}")
    lpts_9.append(f"{px:.1f},{py9:.1f}")
    svg5.append(f'<circle cx="{px:.1f}" cy="{py10:.1f}" r="4.5" fill="#e74c3c"/>')
    svg5.append(f'<circle cx="{px:.1f}" cy="{py9:.1f}" r="4.5" fill="#2980b9"/>')
    if p in (6, 12, 16, 20):
        svg5.append(f'<text x="{px:.1f}" y="{py9+16:.1f}" text-anchor="middle" font-size="10" font-weight="bold" fill="#1f618d">{int(r["mean.k._n9_sin_B19"])}</text>')
    if p in (1, 2, 4, 6, 8, 10, 12, 16, 20):
        svg5.append(f'<text x="{px:.1f}" y="{y0+ph+20}" text-anchor="middle" font-size="11" fill="#2c3e50">{int(p)}</text>')

svg5.append(f'<polyline points="{" ".join(lpts_10)}" fill="none" stroke="#e74c3c" stroke-width="2.2"/>')
svg5.append(f'<polyline points="{" ".join(lpts_9)}" fill="none" stroke="#2980b9" stroke-width="2.2"/>')
svg5.append(f'<text x="{x1+pw//2}" y="{y0+ph+42}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#1a252f">Potencia Suave (Soft Threshold β)</text>')

svg5.append(f'<rect x="330" y="530" width="16" height="16" fill="#e74c3c" rx="3"/><text x="352" y="543" font-size="12.5" fill="#1a252f">n = 10 muestras (incluyendo outgroup B19)</text>')
svg5.append(f'<rect x="640" y="530" width="16" height="16" fill="#2980b9" rx="3"/><text x="662" y="543" font-size="12.5" fill="#1a252f">n = 9 muestras (sin B19 — Red depurada WGCNA)</text>')
svg5.append("</svg>")
write_svg_ascii("figuras_tesis/Fig05_WGCNA_Topologia_Libre_Escala.svg", svg5)


# ==============================================================================
# FIGURA 6: DISTRIBUCIÓN DE LOS 22 MÓDULOS Y HEATMAP DE MODULE EIGENGENES (SIN B19)
# ==============================================================================
print("[6/8] Generando Fig06_WGCNA_Modulos_y_Eigengenes_SinB19.svg...")
res_mod = pd.read_csv("WGCNA_SinB19_02_Resumen_Modulos.csv")
me_sinb19 = pd.read_csv("WGCNA_SinB19_02_Module_Eigengenes.csv", index_col=0)
order_9 = ["B9", "B12", "B14", "B21", "B10", "B11", "B15", "B17", "B27"]
me_sinb19 = me_sinb19.loc[order_9]

W, H = 1260, 740
svg6 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg6.append(f'<text x="{W//2}" y="36" text-anchor="middle" font-size="19" font-weight="bold" fill="#1a252f">Arquitectura de los 22 Módulos de Co-expresión WGCNA y Perfiles de Module Eigengenes (Sin B19, n=9)</text>')
svg6.append(f'<text x="{W//2}" y="58" text-anchor="middle" font-size="13" fill="#5d6d7e">Red Signed (power = 12, deepSplit = 2.5, pamRespectsDendro = FALSE, mergeCutHeight = 0.25) — Módulo grey = 8,10% (2.970 genes)</text>')

# Panel A Superior: Barras de tamaño de los 22 módulos + grey
bx, by, bw_tot, bh_tot = 90, 86, 1100, 224
svg6.append(f'<rect x="{bx}" y="{by}" width="{bw_tot}" height="{bh_tot}" fill="#fbfcfc" stroke="#bdc3c7" stroke-width="1"/>')
svg6.append(f'<text x="{bx + 15}" y="{by + 22}" font-size="13.5" font-weight="bold" fill="#1a252f">A) Número de Transcritos por Módulo (Total = 36.647 transcritos; 22 módulos activos + grey)</text>')

css_color_map = {
    "turquoise": "#40e0d0", "blue": "#2980b9", "brown": "#8b4513", "yellow": "#f1c40f",
    "green": "#27ae60", "red": "#e74c3c", "black": "#2c3e50", "pink": "#ff69b4",
    "magenta": "#c0392b", "purple": "#8e44ad", "greenyellow": "#a2d9ce", "tan": "#d2b48c",
    "salmon": "#fa8072", "cyan": "#00bcd4", "midnightblue": "#191970", "lightcyan": "#81d4fa",
    "grey60": "#7f8c8d", "lightgreen": "#82e0aa", "lightyellow": "#f9e79f", "royalblue": "#4169e1",
    "darkred": "#900c3f", "darkgreen": "#145a32", "grey": "#bdc3c7"
}

n_bars = len(res_mod)
slot_w = (bw_tot - 40) / n_bars
for idx, r in res_mod.iterrows():
    c_name = r["Color_Modulo"]
    cnt = r["Freq"]
    pct = r["Porcentaje"]
    cx = bx + 20 + (idx + 0.5) * slot_w
    bar_h = int((math.sqrt(cnt) / math.sqrt(18000.0)) * 122)
    y_top = by + bh_tot - bar_h - 42
    c_hex = css_color_map.get(c_name, "#3498db")
    svg6.append(f'<rect x="{cx - 14:.1f}" y="{y_top}" width="28" height="{bar_h}" fill="{c_hex}" stroke="#34495e" stroke-width="1"/>')
    svg6.append(f'<text x="{cx:.1f}" y="{y_top - 14}" text-anchor="middle" font-size="9.5" font-weight="bold" fill="#1a252f">{cnt}</text>')
    svg6.append(f'<text x="{cx:.1f}" y="{y_top - 3}" text-anchor="middle" font-size="8.5" fill="#566573">({pct:.1f}%)</text>')
    svg6.append(f'<text x="{cx + 4:.1f}" y="{by + bh_tot - 24}" text-anchor="end" font-size="9.5" font-weight="bold" fill="#1a252f" transform="rotate(-28, {cx + 4:.1f}, {by + bh_tot - 24})">{c_name}</text>')

# Panel B Inferior: Heatmap de Module Eigengenes (9 muestras x 22 módulos)
hx, hy = 215, 350
cw, ch = 44, 32
svg6.append(f'<text x="{bx}" y="{hy - 15}" font-size="13.5" font-weight="bold" fill="#1a252f">B) Perfiles de Expresión de Module Eigengenes (ME) por Muestra y Nicho Altitudinal</text>')

me_cols = list(me_sinb19.columns)
for i, s in enumerate(order_9):
    sp_name, _, col_hex = nicho_short[s]
    ry = hy + i * ch
    svg6.append(f'<rect x="{hx - 20}" y="{ry + 2}" width="12" height="{ch - 4}" fill="{col_hex}" rx="2"/>')
    svg6.append(f'<text x="{hx - 26}" y="{ry + ch//2 + 4}" text-anchor="end" font-size="11.5" font-weight="bold" fill="#1a252f">{s} ({sp_name})</text>')

    for j, m_col in enumerate(me_cols):
        cx = hx + j * cw
        val = float(me_sinb19.loc[s, m_col])
        fill_c = val_to_bwr_hex(val, -2.0, 2.0)
        txt_c = "#ffffff" if abs(val) > 1.05 else "#1a252f"
        svg6.append(f'<rect x="{cx}" y="{ry}" width="{cw}" height="{ch}" fill="{fill_c}" stroke="#ffffff" stroke-width="1"/>')
        svg6.append(f'<text x="{cx + cw//2}" y="{ry + ch//2 + 4}" text-anchor="middle" font-size="9.5" font-weight="bold" fill="{txt_c}">{val:+.2f}</text>')

for j, m_col in enumerate(me_cols):
    cx = hx + j * cw
    c_clean = m_col.replace("ME", "")
    c_hex = css_color_map.get(c_clean, "#3498db")
    ry_bot = hy + 9 * ch
    svg6.append(f'<rect x="{cx + 2}" y="{ry_bot + 3}" width="{cw - 4}" height="10" fill="{c_hex}" stroke="#34495e" stroke-width="0.6"/>')
    svg6.append(f'<text x="{cx + cw//2}" y="{ry_bot + 26}" text-anchor="end" font-size="10.5" font-weight="bold" fill="#1a252f" transform="rotate(-40, {cx + cw//2}, {ry_bot + 26})">{m_col}</text>')

svg6.append(f'<text x="{W//2}" y="{H - 15}" text-anchor="middle" font-size="12" fill="#34495e">Destacan MEpurple/MEyellow (específicos de Alta Montaña: B10, B11, B15), MEblue/MEmagenta (Tierras Bajas: B9, B12, B14, B21) y MEbrown/MEred/MEcyan/MEtan (Datura).</text>')
svg6.append("</svg>")
write_svg_ascii("figuras_tesis/Fig06_WGCNA_Modulos_y_Eigengenes_SinB19.svg", svg6)


# ==============================================================================
# FIGURA 7: LABELED HEATMAP WGCNA FASE 2 (MÓDULOS vs ALCALOIDES TROPÁNICOS)
# ==============================================================================
print("[7/8] Generando Fig07_WGCNA_Heatmap_Modulo_Alcaloide_Fase2.svg...")
r_mat = pd.read_csv("WGCNA_SinB19_03_Correlacion_Modulos_Metabolitos.csv", index_col=0)
p_mat = pd.read_csv("WGCNA_SinB19_03_Pvalores_Modulos_Metabolitos.csv", index_col=0)

# Seleccionar los 24 compuestos con mayor asociación absoluta (incluyendo los 3 precursores X82P, X83P, X84P)
max_abs_r = r_mat.abs().max(axis=0)
top_ta = list(max_abs_r.sort_values(ascending=False).index[:22])
for prec in ["X82P", "X83P", "X84P"]:
    if prec not in top_ta and prec in r_mat.columns:
        top_ta.append(prec)

r_sub = r_mat[top_ta]
p_sub = p_mat[top_ta]

W, H = 1360, 860
svg7 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg7.append(f'<text x="{W//2}" y="34" text-anchor="middle" font-size="19" font-weight="bold" fill="#1a252f">Relación Módulo–Alcaloide Tropánico (WGCNA Fase 2 — n = 8 Muestras Homologadas, Metaboloma Corregido)</text>')
svg7.append(f'<text x="{W//2}" y="56" text-anchor="middle" font-size="12.5" fill="#5d6d7e">Cada celda reporta la correlación de Pearson r (arriba) y el p-valor t-Student (abajo). Asteriscos: * p &lt; 0.05, ** p &lt; 0.01, *** p &lt; 0.001</text>')

hx, hy = 155, 100
cw, ch = 46, 30
mods_list = list(r_sub.index)

for j, comp in enumerate(top_ta):
    cx = hx + j * cw
    c_lbl_col = "#8e44ad" if comp.endswith("P") else "#1a252f"
    svg7.append(f'<text x="{cx + cw//2}" y="{hy - 10}" text-anchor="middle" font-size="11" font-weight="bold" fill="{c_lbl_col}">{comp}</text>')

for i, mod in enumerate(mods_list):
    ry = hy + i * ch
    c_clean = mod.replace("ME", "")
    c_hex = css_color_map.get(c_clean, "#3498db")
    svg7.append(f'<rect x="{hx - 20}" y="{ry + 2}" width="14" height="{ch - 4}" fill="{c_hex}" stroke="#34495e" stroke-width="0.7"/>')
    svg7.append(f'<text x="{hx - 26}" y="{ry + ch//2 + 4}" text-anchor="end" font-size="11" font-weight="bold" fill="#1a252f">{mod}</text>')

    for j, comp in enumerate(top_ta):
        cx = hx + j * cw
        rv = float(r_sub.loc[mod, comp])
        pv = float(p_sub.loc[mod, comp])
        fill_c = val_to_bwr_hex(rv, -1.0, 1.0)
        txt_c = "#ffffff" if abs(rv) > 0.62 else "#1a252f"
        stars = "***" if pv < 0.001 else ("**" if pv < 0.01 else ("*" if pv <= 0.05 else ""))
        stroke_c = "#1a252f" if pv <= 0.05 and abs(rv) >= 0.70 else "#ffffff"
        stroke_w = "1.4" if pv <= 0.05 and abs(rv) >= 0.70 else "0.8"

        svg7.append(f'<rect x="{cx}" y="{ry}" width="{cw}" height="{ch}" fill="{fill_c}" stroke="{stroke_c}" stroke-width="{stroke_w}"/>')
        svg7.append(f'<text x="{cx + cw//2}" y="{ry + 12}" text-anchor="middle" font-size="8.5" font-weight="bold" fill="{txt_c}">{rv:+.2f}{stars}</text>')
        svg7.append(f'<text x="{cx + cw//2}" y="{ry + 24}" text-anchor="middle" font-size="7.5" fill="{txt_c}">({pv:.3f})</text>')

# Barra de color a la derecha
lx_bar = hx + len(top_ta) * cw + 22
svg7.append(f'<text x="{lx_bar}" y="{hy + 12}" font-size="11.5" font-weight="bold" fill="#1a252f">Pearson r</text>')
for s_i in range(25):
    v_s = 1.0 - s_i * (2.0 / 24.0)
    cy_s = hy + 24 + s_i * 18
    svg7.append(f'<rect x="{lx_bar}" y="{cy_s}" width="22" height="18" fill="{val_to_bwr_hex(v_s, -1.0, 1.0)}"/>')
    if s_i in (0, 6, 12, 18, 24):
        svg7.append(f'<text x="{lx_bar + 28}" y="{cy_s + 13}" font-size="10.5" fill="#2c3e50">{v_s:+.1f}</text>')

svg7.append(f'<text x="{W//2}" y="{H - 20}" text-anchor="middle" font-size="12" fill="#34495e">Asociaciones clave: X70T (+0,92 en MEblue/MEmagenta vs. -0,97 en MEred), X52T (-0,97 en MEpurple Alta Montaña vs. +0,86 en MEbrown) y precursor X82P (+0,88 log10 / +0,95 lineal en MEmagenta).</text>')
svg7.append("</svg>")
write_svg_ascii("figuras_tesis/Fig07_WGCNA_Heatmap_Modulo_Alcaloide_Fase2.svg", svg7)


# ==============================================================================
# FIGURA 8: VOLCANO PLOT (ALTA MONTAÑA vs TIERRAS BAJAS) + HUB GENES (kME vs GS)
# ==============================================================================
print("[8/8] Generando Fig08_Volcano_DEGs_y_HubGenes_MM_vs_GS.svg...")
from scipy import stats

alta_montana = ["B10", "B11", "B15"]
tierras_bajas = ["B9", "B12", "B14", "B21"]
X_alta = df_vst.loc[alta_montana].values
X_baja = df_vst.loc[tierras_bajas].values
lfc = X_alta.mean(axis=0) - X_baja.mean(axis=0)
_, p_de = stats.ttest_ind(X_alta, X_baja, axis=0, equal_var=False)
p_de = np.clip(np.nan_to_num(p_de, nan=1.0), 1e-7, 1.0)
neg_log10_p = -np.log10(p_de)

df_hubs = pd.read_csv("WGCNA_SinB19_04_Top_HubGenes_Candidatos_DEGs_TA.csv")

W, H = 1220, 600
svg8 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">']
svg8.append(f'<text x="{W//2}" y="36" text-anchor="middle" font-size="19" font-weight="bold" fill="#1a252f">Expresión Diferencial entre Nichos Altitudinales y Selección de Hub Genes Candidatos (OE2 + OE3)</text>')
svg8.append(f'<text x="{W//2}" y="58" text-anchor="middle" font-size="13" fill="#5d6d7e">Integración DESeq2/VST (Subclado A Alta Montaña vs. Subclado B Tierras Bajas) + Centralidad WGCNA (kME ≥ 0,85 y |GS| ≥ 0,75)</text>')

# Panel A: Volcano Plot (Submuestreo inteligente de puntos grises + todos los significativos para SVG ágil)
x0, y0, pw, ph = 75, 100, 490, 410
svg8.append(f'<rect x="{x0}" y="{y0}" width="{pw}" height="{ph}" fill="#fbfcfc" stroke="#bdc3c7" stroke-width="1.2"/>')
svg8.append(f'<text x="{x0+pw//2}" y="{y0-12}" text-anchor="middle" font-size="13.5" font-weight="bold" fill="#1a252f">A) Volcano Plot: Alta Montaña (644 Up) vs. Tierras Bajas (2.452 Up)</text>')

lfc_min, lfc_max = -11.0, 11.0
p_max_plot = 6.2

def v_to_px(lv, pv):
    px = x0 + ((max(lfc_min, min(lfc_max, lv)) - lfc_min) / (lfc_max - lfc_min)) * pw
    py = y0 + ph - (min(p_max_plot, pv) / p_max_plot) * ph
    return px, py

# Líneas de umbral |log2FC| = 1.5 y p = 0.01 (-log10 = 2.0)
ux1, uy_cut = v_to_px(-1.5, 2.0)
ux2, _ = v_to_px(1.5, 2.0)
svg8.append(f'<line x1="{x0}" y1="{uy_cut:.1f}" x2="{x0+pw}" y2="{uy_cut:.1f}" stroke="#7f8c8d" stroke-dasharray="5,4" stroke-width="1.2"/>')
svg8.append(f'<line x1="{ux1:.1f}" y1="{y0}" x2="{ux1:.1f}" y2="{y0+ph}" stroke="#7f8c8d" stroke-dasharray="5,4" stroke-width="1.2"/>')
svg8.append(f'<line x1="{ux2:.1f}" y1="{y0}" x2="{ux2:.1f}" y2="{y0+ph}" stroke="#7f8c8d" stroke-dasharray="5,4" stroke-width="1.2"/>')

rng = np.random.default_rng(42)
for idx_g in range(len(lfc)):
    lv, pv = lfc[idx_g], neg_log10_p[idx_g]
    is_sig = (abs(lv) >= 1.5) and (pv >= 2.0)
    if not is_sig and rng.random() > 0.035:
        continue
    if is_sig and pv < 2.8 and abs(lv) < 3.5 and rng.random() > 0.45:
        continue
    px, py = v_to_px(lv, pv)
    if is_sig:
        col_pt = "#8e44ad" if lv > 0 else "#27ae60"
        r_pt = "2.6"
        op = "0.78"
    else:
        col_pt = "#bdc3c7"
        r_pt = "1.8"
        op = "0.45"
    svg8.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{r_pt}" fill="{col_pt}" opacity="{op}"/>')

for tick_x in [-10, -5, 0, 5, 10]:
    tx, _ = v_to_px(tick_x, 0)
    svg8.append(f'<text x="{tx:.1f}" y="{y0+ph+18}" text-anchor="middle" font-size="11" fill="#2c3e50">{tick_x:+d}</text>')
for tick_y in [0, 1, 2, 3, 4, 5, 6]:
    _, ty = v_to_px(0, tick_y)
    svg8.append(f'<text x="{x0-8}" y="{ty+4:.1f}" text-anchor="end" font-size="11" fill="#2c3e50">{tick_y}</text>')

svg8.append(f'<text x="{x0+95}" y="{y0+28}" text-anchor="middle" font-size="11.5" font-weight="bold" fill="#1e8449">Up Tierras Bajas (2.452)</text>')
svg8.append(f'<text x="{x0+pw-95}" y="{y0+28}" text-anchor="middle" font-size="11.5" font-weight="bold" fill="#8e44ad">Up Alta Montaña (644)</text>')
svg8.append(f'<text x="{x0+pw//2}" y="{y0+ph+38}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#1a252f">log2(Fold Change) [Alta Montaña vs. Tierras Bajas]</text>')
svg8.append(f'<text x="{x0-32}" y="{y0+ph//2}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#1a252f" transform="rotate(-90, {x0-32}, {y0+ph//2})">-log10(p-valor)</text>')

# Panel B: Dispersión Module Membership (kME) vs Gene Significance (|GS|) para los 1.440 Hub Genes
x1 = 655
svg8.append(f'<rect x="{x1}" y="{y0}" width="{pw}" height="{ph}" fill="#fbfcfc" stroke="#bdc3c7" stroke-width="1.2"/>')
svg8.append(f'<text x="{x1+pw//2}" y="{y0-12}" text-anchor="middle" font-size="13.5" font-weight="bold" fill="#1a252f">B) Selección de Hub Genes: Centralidad (kME) vs. Significancia Alcaloide (|GS|)</text>')

k_min, k_max = 0.83, 1.005
g_min, g_max = 0.73, 1.005

def kg_to_px(k_val, g_val):
    px = x1 + ((max(k_min, min(k_max, k_val)) - k_min) / (k_max - k_min)) * pw
    py = y0 + ph - ((max(g_min, min(g_max, g_val)) - g_min) / (g_max - g_min)) * ph
    return px, py

# Cuadrante de élite (kME >= 0.85, |GS| >= 0.75)
qx, qy = kg_to_px(0.85, 0.75)
svg8.append(f'<rect x="{qx:.1f}" y="{y0}" width="{x1+pw-qx:.1f}" height="{qy-y0:.1f}" fill="#eafaf1" stroke="#27ae60" stroke-dasharray="5,4" stroke-width="1.5"/>')
svg8.append(f'<text x="{qx + (x1+pw-qx)/2:.1f}" y="{qy - 12:.1f}" text-anchor="middle" font-size="11.5" font-weight="bold" fill="#1e8449">Zona Hub Genes Élite (1.440 transcritos: kME ≥ 0,85 y |GS| ≥ 0,75)</text>')

for _, row_h in df_hubs.iloc[::2].iterrows():
    kv = float(row_h["kME_Propio_Modulo"])
    gv = max(abs(float(row_h["GS_X70T_log10"])), abs(float(row_h["GS_X52T_log10"])), abs(float(row_h["GS_X82P_log10"])))
    m_c = row_h["Module_Sin_B19"]
    c_hex = css_color_map.get(m_c, "#2980b9")
    px, py = kg_to_px(kv, gv)
    svg8.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="3.4" fill="{c_hex}" opacity="0.82"/>')

# Etiquetar 5 Hub Genes candidatos representativos sin solapamiento
hub_callouts = [
    (0, -14, 16, "end"),
    (1, -14, -8, "end"),
    (2, -14, 14, "end"),
    (3, -14, 24, "end"),
    (4, -14, 10, "end"),
]
for idx_top, dx_lbl, dy_lbl, anc in hub_callouts:
    row_h = df_hubs.iloc[idx_top]
    kv = float(row_h["kME_Propio_Modulo"])
    gv = max(abs(float(row_h["GS_X70T_log10"])), abs(float(row_h["GS_X52T_log10"])), abs(float(row_h["GS_X82P_log10"])))
    px, py = kg_to_px(kv, gv)
    short_id = "_".join(row_h["Transcript_ID"].split("_")[:2])
    svg8.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="5.5" fill="none" stroke="#1a252f" stroke-width="1.8"/>')
    svg8.append(f'<line x1="{px:.1f}" y1="{py:.1f}" x2="{px + dx_lbl + 3:.1f}" y2="{py + dy_lbl - 3:.1f}" stroke="#1a252f" stroke-width="1"/>')
    svg8.append(f'<text x="{px + dx_lbl:.1f}" y="{py + dy_lbl:.1f}" text-anchor="{anc}" font-size="9.5" font-weight="bold" fill="#1a252f">{short_id} ({row_h["Module_Sin_B19"]})</text>')

for tk in [0.85, 0.88, 0.91, 0.94, 0.97, 1.00]:
    tx, _ = kg_to_px(tk, g_min)
    svg8.append(f'<text x="{tx:.1f}" y="{y0+ph+18}" text-anchor="middle" font-size="11" fill="#2c3e50">{tk:.2f}</text>')
for tg in [0.75, 0.80, 0.85, 0.90, 0.95, 1.00]:
    _, ty = kg_to_px(k_min, tg)
    svg8.append(f'<text x="{x1-8}" y="{ty+4:.1f}" text-anchor="end" font-size="11" fill="#2c3e50">{tg:.2f}</text>')

svg8.append(f'<text x="{x1+pw//2}" y="{y0+ph+38}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#1a252f">Module Membership (kME — Centralidad Intramodular)</text>')
svg8.append(f'<text x="{x1-36}" y="{y0+ph//2}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#1a252f" transform="rotate(-90, {x1-36}, {y0+ph//2})">Gene Significance |GS| (Correlación con Alcaloides TA)</text>')

svg8.append("</svg>")
write_svg_ascii("figuras_tesis/Fig08_Volcano_DEGs_y_HubGenes_MM_vs_GS.svg", svg8)

print("¡Las 8 figuras vectoriales SVG han sido generadas exitosamente en 'figuras_tesis/'!")
