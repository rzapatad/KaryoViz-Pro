import streamlit as st
from pycirclize import Circos
import pandas as pd
import matplotlib.pyplot as plt
import io, os, re, random

# --- 1. CONFIGURACIÓN DE PÁGINA ---
st.set_page_config(page_title="KaryoViz Pro - Richard Zapata", layout="wide", initial_sidebar_state="expanded")

# --- 2. ESTILOS CSS ---
st.markdown("""
    <style>
        .block-container { padding-top: 1rem; padding-bottom: 0rem; }
        section[data-testid="stSidebar"] [data-testid="stTable"],
        section[data-testid="stSidebar"] [data-testid="stDataFrame"] { width: 100% !important; }
        #MainMenu, footer, header {visibility: hidden;}
        .pro-title {text-align: center; font-family: 'Inter', sans-serif; font-weight: 800; color: #1E1E1E; letter-spacing: -1px; margin-bottom: -5px;}
        .pro-subtitle {text-align: center; font-family: 'Inter', sans-serif; color: #666; font-size: 0.95em; padding-bottom: 20px;}
        .stButton>button {border-radius: 8px; font-weight: 600; transition: all 0.3s;}
        .stMetric { background-color: #f8f9fa; padding: 10px; border-radius: 10px; border: 1px solid #eee; text-align: center;}
    </style>
    <h1 class="pro-title">🧬 KaryoViz Pro</h1>
    <p class="pro-subtitle">Advanced Cytogenetic Mapping & Structural Visualization</p>
""", unsafe_allow_html=True)

# --- 3. GESTIÓN DE ESTADO ---
def cargar_demo():
    st.session_state["text_l_manual"] = (
        "chr1,start1,end1,chr2,start2,end2,color,fusion_gene\n"
        "chr9,133500000,133700000,chr22,23200000,23400000,#E74C3C,BCR::ABL1\n"
        "chr15,20000000,21000000,chr17,38000000,39000000,#3498DB,PML::RARA\n"
        "chr8,127735433,127742951,chr14,105586437,106032614,#F1C40F,MYC::IGH\n"
        "chr12,12000000,13000000,chr21,35000000,36000000,#9B59B6,ETV6::RUNX1\n"
        "chr11,118000000,119000000,chr4,153000000,154000000,#2ECC71,KMT2A::AFF1\n"
        "chr2,29000000,30000000,chr5,170000000,171000000,#E67E22,ALK::EML4\n"
        "chr14,106000000,107000000,chr18,60000000,61000000,#D35400,BCL2::IGH\n"
        "chr6,167000000,168000000,chr9,130000000,131000000,#34495E,DEK::NUP214\n"
        "chr3,52000000,53000000,chr21,44000000,45000000,#7F8C8D,RUNX1::MECOM\n"
        "chr1,153000000,154000000,chr19,53000000,54000000,#1ABC9C,TP53::MDM2"
    )
    st.session_state["text_c_manual"] = (
        "chromosome,start,end,value\nchr7,0,159345973,1.0\nchr13,0,114364328,-1.0\nchr21,0,46709983,1.2\n"
        "chr1,120000000,150000000,-0.8\nchr8,0,145138636,0.9\nchr18,0,80373285,-0.7\nchrX,0,156040895,1.5"
    )
    st.session_state["text_s_manual"] = (
        "chromosome,start,value\nchr1,50000000,40\nchr17,20000000,95\nchr22,23000000,88\nchr2,29445000,100"
    )

def nuevo_analisis():
    for key in ['text_l_manual', 'text_c_manual', 'text_s_manual', 'df_l_final', 'df_c_final', 'iscn_val']:
        st.session_state[key] = "" if 'text' in key or 'iscn' in key else None
    st.session_state.file_uploader_key = st.session_state.get('file_uploader_key', 0) + 1

# Inicialización
for k in ["text_l_manual", "text_c_manual", "text_s_manual", "iscn_val"]:
    if k not in st.session_state: st.session_state[k] = ""
if 'df_l_final' not in st.session_state: st.session_state.df_l_final = None
if 'df_c_final' not in st.session_state: st.session_state.df_c_final = None
if "file_uploader_key" not in st.session_state: st.session_state.file_uploader_key = 0

CANCER_GENES = ["RUNX1", "RUNX1T1", "BCR", "ABL1", "PML", "RARA", "MYC", "IGH", "BCL2", "ALK", "EML4", "KMT2A", "AFF1", "DEK", "NUP210"]

# --- 4. BARRA LATERAL ---
with st.sidebar:
    tabs = st.tabs(["🎨 Diseño & Data", "🔍 Analizador ISCN"])
    with tabs[0]:
        c1, c2, c3 = st.columns(3)
        with c1: st.button("🧹 Nuevo", use_container_width=True, on_click=nuevo_analisis)
        with c2: st.button("🚀 Demo", type="primary", use_container_width=True, on_click=cargar_demo)
        st.divider()
        with st.expander("🖼️ Ajustes Visuales", expanded=True):
            genome_ver = st.selectbox("Referencia", ["hg38", "hg19"])
            zoom_scale = st.slider("Canvas", 5, 20, 12)
            label_size = st.slider("Fuente (Chr)", 6, 16, 10)
            show_hotspots = st.checkbox("Heatmap Hotspots", value=True)
            h_sens = st.slider("Sensibilidad", 1, 10, 5)
        with st.expander("🔗 Links"):
            up_l_man = st.file_uploader("Subir CSV", type="csv", key=f"l_up_{st.session_state.file_uploader_key}")
            st.text_area("Pegar Links", key="text_l_manual", height=80)
            l_col_ui, l_wid = st.color_picker("Color Base", "#1f77b4"), st.slider("Grosor", 0.1, 5.0, 1.2)
            l_opa = st.slider("Opacidad", 0.1, 1.0, 0.6)
            f_cols = st.columns(2)
            with f_cols[0]: tag_size = st.number_input("Tag Size", 5, 20, 9)
            with f_cols[1]: tag_color = st.color_picker("Tag Color", "#D35400")
        with st.expander("📊 CNV"):
            up_c = st.file_uploader("Subir CSV", type="csv", key=f"c_up_{st.session_state.file_uploader_key}")
            st.text_area("Pegar CNV", key="text_c_manual", height=80)
            c_cols = st.columns(2)
            with c_cols[0]: cp = st.color_picker("Ganancia", "#77DD77")
            with c_cols[1]: cn = st.color_picker("Pérdida", "#FF6961")
        with st.expander("📍 SNPs"):
            up_s = st.file_uploader("Subir CSV", type="csv", key=f"s_up_{st.session_state.file_uploader_key}")
            st.text_area("Pegar SNPs", key="text_s_manual", height=80)
            snp_col, snp_width_f = st.color_picker("Color SNP", "#0000FF"), st.slider("Grosor SNP", 50, 500, 200)

    with tabs[1]:
        st.header("Analizador ISCN")
        with st.expander("📋 Demos Clínicas (40 Casos)", expanded=True):
            demo_opts = {
                "Leucemias (10 ej)": "t(9;22)(q34;q11.2)\nt(15;17)(q24.2;q21.2)\nt(8;21)(q22;q22)\ninv(16)(p13.1;q22)\nt(12;21)(p13.2;q22.1)\nt(11;19)(q23.3;p13.3)\nt(1;19)(q23.3;p13.3)\nt(6;9)(p23;q34.1)\nt(9;11)(p22;q23.3)\ninv(3)(q21.3;q26.2)",
                "Linfomas (10 ej)": "t(8;14)(q24.21;q32.33)\nt(14;18)(q32.3;q21.3)\nt(11;14)(q13.3;q32.3)\nt(2;8)(p11.2;q24.21)\nt(8;22)(q24.21;q11.2)\nt(3;14)(q27.3;q32.33)\nt(11;18)(q21.31;q21.1)\nt(14;19)(q32.3;q13.3)\nt(9;14)(p13.2;q32.33)\nt(10;14)(q24.31;q32.3)",
                "Mielomas (10 ej)": "t(4;14)(p16.3;q32.3)\nt(14;16)(q32.3;q23.3)\nt(14;20)(q32.3;q12.2)\ndel(17)(p13)\n+1(q21)\ndel(13)(q14)\nt(6;14)(p21.1;q32.3)\n+9\n+11\n+15",
                "Paneles moleculares (10 ej)": "BCR::ABL1\nPML::RARA\nRUNX1::RUNX1T1\nEML4::ALK\nKMT2A::AFF1\nMYC::IGH\nBCL2::IGH\nPAX5::JAK2\nETV6::NTRK3\nNUP98::NSD1"
            }
            sel_demo = st.selectbox("Seleccionar set:", list(demo_opts.keys()))
            if st.button("🚀 Aplicar Demo"):
                st.session_state.iscn_val = demo_opts[sel_demo]
                st.rerun()
        st.text_area("Cariotipos / Genes", height=150, key="iscn_val")
        up_iscn_p2 = st.file_uploader("Adjuntar lista (CSV)", type="csv", key=f"iscn_up_{st.session_state.file_uploader_key}")
        b1, b2 = st.columns(2)
        with b1: btn_proc = st.button("🚀 Procesar Fórmula", type="primary", use_container_width=True)
        with b2: st.button("Limpiar", use_container_width=True, on_click=nuevo_analisis)
        if st.session_state.df_l_final is not None and not st.session_state.df_l_final.empty:
            st.subheader("🧬 Hallazgos Estructurales")
            st.session_state.df_l_final = st.data_editor(st.session_state.df_l_final, hide_index=True, key="ed_est")
        if st.session_state.df_c_final is not None and not st.session_state.df_c_final.empty:
            st.subheader("📊 Hallazgos Numéricos (CNV)")
            st.session_state.df_c_final = st.data_editor(st.session_state.df_c_final, hide_index=True, key="ed_num")

# --- 5. MOTOR DE PROCESAMIENTO ---
base_p = os.path.dirname(__file__)
size_f = os.path.join(base_p, f"{genome_ver}.chrom.sizes.txt")
band_f = os.path.join(base_p, f"{genome_ver}_cytoBand.txt")
genes_db_f = os.path.join(base_p, f"{genome_ver}_genes.csv")

def find_best_gene(chrom, start, end, df_genes):
    if df_genes is None or df_genes.empty: return ""
    m = df_genes[(df_genes["chrom"] == chrom) & (df_genes["Start"] < end) & (df_genes["End"] > start)]
    for cg in CANCER_GENES:
        if not m.empty and cg in m["geneSymbol"].values: return cg
    return m.iloc[0]["geneSymbol"] if not m.empty else ""

def procesar_dual_logic(lista, df_bandas, df_tamanos, df_genes=None):
    enlaces, numericas = [], []
    df_bandas = df_bandas.dropna(subset=["band"])
    for item in lista:
        if pd.isna(item) or not str(item).strip(): continue
        txt = str(item).replace(" ", "").upper()
        if "::" in txt and df_genes is not None:
            pts = txt.split("::")
            g1, g2 = df_genes[df_genes["geneSymbol"].str.upper() == pts[0]], df_genes[df_genes["geneSymbol"].str.upper() == pts[1]]
            if not g1.empty and not g2.empty:
                enlaces.append({"chr1": g1.iloc[0,0], "start1": int(g1.iloc[0,1]), "end1": int(g1.iloc[0,2]), "chr2": g2.iloc[0,0], "start2": int(g2.iloc[0,1]), "end2": int(g2.iloc[0,2]), "color": "#E74C3C", "fusion_gene": txt})
            continue
        e_m = re.search(r"(T|DEL|DUP|INV|ROB|R|ISO|ADD|INS|[\+\-])\(?(.*?)\)?\((.*?)\)", txt)
        if e_m:
            tipo, c, b = e_m.groups()
            c_l, b_l = c.split(";"), (b.split(";") if ";" in b else re.findall(r'[PQ]\d+', b))
            if tipo in ["DEL", "ADD", "INS", "DUP", "+", "-"]:
                chr_cnv = f"chr{c_l[0]}"
                r_cnv = df_bandas[(df_bandas["chr"] == chr_cnv) & (df_bandas["band"].str.contains(b_l[0], case=False))].iloc[:1]
                if not r_cnv.empty:
                    val = -1.0 if tipo in ["DEL", "-"] else 1.0
                    numericas.append({"chromosome": chr_cnv, "start": r_cnv.iloc[0]["start"], "end": r_cnv.iloc[0]["end"], "value": val})
            elif tipo in ["T", "INV", "ROB", "R", "ISO"]:
                if len(c_l) >= 1 and len(b_l) >= 2:
                    c1, c2 = f"chr{c_l[0]}", (f"chr{c_l[1]}" if len(c_l)>1 else f"chr{c_l[0]}")
                    r1 = df_bandas[(df_bandas["chr"] == c1) & (df_bandas["band"].str.contains(b_l[0], case=False))].iloc[:1]
                    r2 = df_bandas[(df_bandas["chr"] == c2) & (df_bandas["band"].str.contains(b_l[1], case=False))].iloc[:1]
                    if not r1.empty and not r2.empty:
                        gn1, gn2 = find_best_gene(c1, r1.iloc[0]["start"], r1.iloc[0]["end"], df_genes), find_best_gene(c2, r2.iloc[0]["start"], r2.iloc[0]["end"], df_genes)
                        enlaces.append({"chr1": c1, "start1": r1.iloc[0]["start"], "end1": r1.iloc[0]["end"], "chr2": c2, "start2": r2.iloc[0]["start"], "end2": r2.iloc[0]["end"], "color": random.choice(["#5DADE2", "#A569BD"]), "fusion_gene": f"{gn1}::{gn2}" if gn1 else ""})
            continue
        n_m = re.search(r"([\+-])(\d+|X|Y)", txt)
        if n_m:
            s, ch = n_m.groups(); c_n = f"chr{ch.upper()}"; sz = df_tamanos[df_tamanos["chr"] == c_n]
            if not sz.empty: numericas.append({"chromosome": c_n, "start": 0, "end": sz.iloc[0]["size"], "value": 1.0 if s=="+" else -1.0})
    return pd.DataFrame(enlaces), pd.DataFrame(numericas)

if btn_proc:
    if os.path.exists(band_f) and os.path.exists(size_f):
        db_b, db_s = pd.read_csv(band_f, sep="\t", header=None, names=["chr","start","end","band","gie"]), pd.read_csv(size_f, sep="\t", header=None, names=["chr","size"])
        db_g = pd.read_csv(genes_db_f) if os.path.exists(genes_db_f) else None
        if up_iscn_p2:
            df_csv = pd.read_csv(up_iscn_p2)
            lista_iscn = df_csv[df_csv.columns[0]].astype(str).tolist()
        else:
            lista_iscn = st.session_state.iscn_val.split("\n")
        st.session_state.df_l_final, st.session_state.df_c_final = procesar_dual_logic(lista_iscn, db_b, db_s, db_g)
        st.rerun()

# --- 6. RENDERIZADO VISUAL ---
if os.path.exists(size_f) and os.path.exists(band_f):
    df_sz_raw = pd.read_csv(size_f, sep="\t", header=None, names=["chr","size"])
    logical_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    df_sz = df_sz_raw[df_sz_raw["chr"].isin(logical_order)].copy()
    df_sz['idx'] = df_sz['chr'].map({c: i for i, c in enumerate(logical_order)})
    df_sz = df_sz.sort_values('idx').drop('idx', axis=1)
    chr_max_sizes = {r["chr"]: r["size"] for _, r in df_sz.iterrows()}

    l_pts = [st.session_state.df_l_final, (pd.read_csv(io.StringIO(st.session_state.text_l_manual)) if st.session_state.text_l_manual.strip() else None), (pd.read_csv(up_l_man) if up_l_man else None)]
    d_l = pd.concat([p for p in l_pts if p is not None]).drop_duplicates() if any(p is not None for p in l_pts) else pd.DataFrame()
    c_pts = [st.session_state.df_c_final, (pd.read_csv(io.StringIO(st.session_state.text_c_manual)) if st.session_state.text_c_manual.strip() else None), (pd.read_csv(up_c) if up_c else None)]
    d_c = pd.concat([p for p in c_pts if p is not None]).drop_duplicates() if any(p is not None for p in c_pts) else pd.DataFrame()
    s_pts = [(pd.read_csv(io.StringIO(st.session_state.text_s_manual)) if st.session_state.text_s_manual.strip() else None), (pd.read_csv(up_s) if up_s else None)]
    d_s = pd.concat([p for p in s_pts if p is not None]).drop_duplicates() if any(p is not None for p in s_pts) else pd.DataFrame()

    m1, m2, m3, m4, m5 = st.columns(5)
    with m1: st.metric("Fusiones", len(d_l))
    with m2: 
        h_val = d_l.iloc[:,0].mode()[0].replace("chr","") if not d_l.empty else "-"
        st.metric("Hotspot", f"Chr {h_val}")
    with m3: st.metric("Gain", len(d_c[d_c.iloc[:,3] > 0]) if not d_c.empty else 0)
    with m4: st.metric("Loss", len(d_c[d_c.iloc[:,3] < 0]) if not d_c.empty else 0)
    with m5: st.metric("SNPs", len(d_s))

    st.write("🛠️ **Capas & Búsqueda Genómica**:")
    lay1, lay2, lay3, lay4 = st.columns(4)
    with lay1: show_links = st.toggle("Links", value=True)
    with lay2: show_cnv = st.toggle("CNV", value=True)
    with lay3: show_snp = st.toggle("SNPs", value=True)
    with lay4: show_tags = st.toggle("Tags", value=True)

    gene_search = st.text_input("🔍 Localizar Gen (Ej: ALK, TP53, BCR):", "").upper().strip()

    circos = Circos({r["chr"]: r["size"] for _, r in df_sz.iterrows()}, space=3)
    db_band_raw = pd.read_csv(band_f, sep="\t", header=None, names=["chr","start","end","band","gieStain"])
    gie_c = {"gneg":"#FFFFFF","gpos25":"#D9D9D9","gpos50":"#969696","gpos75":"#525252","gpos100":"#000000","acen":"#800000","stalk":"#707070"}

    df_genes_db = pd.read_csv(genes_db_f) if os.path.exists(genes_db_f) else None

    for sector in circos.sectors:
        max_len = chr_max_sizes.get(sector.name, 0)
        sec_b = db_band_raw[db_band_raw["chr"] == sector.name]
        track_b = sector.add_track((65, 70))
        for _, r in sec_b.iterrows(): track_b.rect(r["start"], min(r["end"], max_len), color=gie_c.get(r["gieStain"], "#FFF"), ec="#DDD", lw=0.3)
        hit_count = len(d_l[d_l.iloc[:,0] == sector.name]) + len(d_l[d_l.iloc[:,3] == sector.name]) if not d_l.empty else 0
        if show_hotspots: sector.add_track((71, 75)).rect(0, sector.size, color=plt.cm.YlOrRd(min(hit_count/h_sens, 1.0)), ec="none", alpha=0.7)
        sector.text(sector.name.replace("chr", ""), r=92, size=label_size, fontweight='bold', color="#333")
        
        if gene_search and df_genes_db is not None:
            found_gene = df_genes_db[(df_genes_db["geneSymbol"].str.upper() == gene_search) & (df_genes_db["chrom"] == sector.name)]
            if not found_gene.empty:
                g_pos = (found_gene.iloc[0]["Start"] + found_gene.iloc[0]["End"]) / 2
                sector.add_track((77, 85)).rect(found_gene.iloc[0]["Start"], found_gene.iloc[0]["End"], color="red", ec="red", lw=2)
                sector.text(gene_search, x=g_pos, r=105, size=label_size+2, color="red", fontweight='bold')

        if show_cnv and not d_c.empty:
            d = d_c[d_c.iloc[:,0] == sector.name]
            if not d.empty:
                track_c = sector.add_track((58, 63))
                for _, rc in d.iterrows():
                    s_pos, e_pos = int(rc.iloc[1]), min(int(rc.iloc[2]), max_len)
                    if s_pos < max_len: track_c.rect(s_pos, e_pos, color=cp if rc.iloc[3] > 0 else cn, alpha=0.9, ec="white", lw=0.5)
        
        if show_snp and not d_s.empty:
            s = d_s[d_s.iloc[:,0] == sector.name]
            if not s.empty:
                track_s = sector.add_track((77, 85))
                for _, rs in s.iterrows():
                    if int(rs.iloc[1]) < max_len: track_s.bar([int(rs.iloc[1])], [rs.iloc[2]], width=sector.size/snp_width_f, color=snp_col, alpha=0.8)
        
        if show_tags and not d_l.empty:
            for _, row in d_l.iterrows():
                if row.iloc[0] == sector.name:
                    l_t = str(row['fusion_gene']) if 'fusion_gene' in d_l.columns else str(row.iloc[-1])
                    if l_t and l_t != "nan": sector.text(l_t.split("::")[0], r=52, x=min(int(row.iloc[1]), max_len), size=tag_size, color=tag_color, fontweight='bold')
    
    if show_links and not d_l.empty:
        for _, r in d_l.iterrows():
            try:
                c1, s1, e1, c2, s2, e2 = r.iloc[0], r.iloc[1], r.iloc[2], r.iloc[3], r.iloc[4], r.iloc[5]
                s1, e1 = min(s1, chr_max_sizes.get(c1, s1)), min(e1, chr_max_sizes.get(c1, e1))
                s2, e2 = min(s2, chr_max_sizes.get(c2, s2)), min(e2, chr_max_sizes.get(c2, e2))
                circos.link((c1, s1, e1), (c2, s2, e2), color=r.get('color', l_col_ui), lw=l_wid, alpha=l_opa)
            except: continue

    fig, ax = plt.subplots(figsize=(zoom_scale, zoom_scale), subplot_kw={'projection': 'polar'})
    circos.plotfig(ax=ax); ax.set_axis_off(); plt.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.07)
    st.pyplot(fig, use_container_width=True)
    
    st.markdown("### 📚 Leyenda Citogenética")
    col_leg1, col_leg2, col_leg3 = st.columns(3)
    with col_leg1:
        st.write("**Bandas (GIE Stains)**")
        st.write("⚪ *gneg*: Claras"); st.write("🔘 *gpos*: Grises"); st.write("⚫ *gpos100*: Oscuras")
    with col_leg2:
        st.write("**Regiones**")
        st.write("🔴 *acen*: Centrómero"); st.write("〰️ *stalk*: Satélites"); st.write("🔥 *Heatmap*: Hotspots")
    with col_leg3:
        st.write("**Marcas**")
        st.write(f"🟩 *Gain*: {cp}"); st.write(f"🟥 *Loss*: {cn}"); st.write(f"🔵 *Links*: Fusiones")

    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=300, transparent=True, bbox_inches='tight', pad_inches=0.4)
    st.download_button(label="📥 Guardar HQ PNG", data=buf.getvalue(), file_name="KaryoViz_HQ.png", mime="image/png", use_container_width=True)
else:
    st.error(f"Archivos de referencia faltantes.")
