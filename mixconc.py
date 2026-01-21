import streamlit as st
import pandas as pd
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Frame, PageTemplate
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
import tempfile
import time

# =========================
# 1. 全局配置与样式
# =========================
APP_VERSION = "v1.2"
st.set_page_config(page_title="液体混合计算器", page_icon="🧪", layout="wide")

# CSS 样式优化
st.markdown("""
    <style>
    /* 0. 顶部留白调整：减少主容器顶部的 padding */
    .block-container {
        padding-top: 2rem !important; /* 默认通常是 6rem 左右，设为 2rem 即可减少约 2/3 */
        padding-bottom: 2rem !important;
    }

    /* 1. 侧边栏宽度调整 (约1.1-1.2倍) */
    [data-testid="stSidebar"] {
        min-width: 400px !important;
        max-width: 400px !important;
    }

    /* 2. 结果卡片样式 */
    .metric-card {
        background-color: #f0f2f6;
        padding: 20px;
        border-radius: 10px;
        border-left: 5px solid #4e8cff;
        box-shadow: 2px 2px 5px rgba(0,0,0,0.05);
        text-align: center;
        margin-bottom: 20px;
    }
    .metric-value {
        font-size: 1.8em;
        font-weight: bold;
        color: #31333F;
    }
    .metric-label {
        font-size: 0.9em;
        color: #666;
        margin-bottom: 5px;
    }
    
    /* 3. 隐藏数字输入框微调按钮 */
    input[type=number]::-webkit-inner-spin-button, 
    input[type=number]::-webkit-outer-spin-button { 
        -webkit-appearance: none; 
        margin: 0; 
    }
    
    /* 4. 参考信息样式 */
    .ref-text {
        font-size: 0.85em;
        color: #555;
        background-color: #eef;
        padding: 8px;
        border-radius: 5px;
        margin-bottom: 15px;
        line-height: 1.6;
    }
    
    /* 5. 强制表格内容居中 (针对 Streamlit 原生表格组件) */
    th {
        text-align: center !important;
    }
    td {
        text-align: center !important;
    }
    </style>
    """, unsafe_allow_html=True)

# =========================
# 2. PDF 字体注册
# =========================
try:
    pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
    FONT_NAME = 'STSong-Light'
except Exception as e:
    FONT_NAME = 'Helvetica'

# =========================
# 3. 核心工具函数
# =========================
MASS_UNIT_TO_G = {"μg": 1e-6, "mg": 1e-3, "g": 1.0, "kg": 1e3}
VOL_UNIT_TO_ML = {"μL": 1e-3, "mL": 1.0, "L": 1000.0}
CONC_MASS_UNIT_TO_G_PER_L = {"μg/L": 1e-6, "mg/L": 1e-3, "g/L": 1.0}

def get_water_density(t):
    return 1.0 - 0.0003 * (t - 4)

def get_saline_density(t):
    return 1.004 - 0.0003 * (t - 20)

def auto_format_solute(mass_g):
    """根据溶质质量大小自动选择单位"""
    if mass_g == 0: return "0.00 g"
    if mass_g < 1e-6: return f"{mass_g * 1e9:.3f} ng"
    elif mass_g < 1e-3: return f"{mass_g * 1e6:.3f} μg"
    elif mass_g < 1.0: return f"{mass_g * 1e3:.3f} mg"
    else: return f"{mass_g:.3f} g"

def calculate_solute_mass(conc, unit, molar_mass, density, total_mass_g):
    """计算溶质绝对质量(g)"""
    if unit in CONC_MASS_UNIT_TO_G_PER_L:
        vol_L = (total_mass_g / density) / 1000.0
        conc_g_L = conc * CONC_MASS_UNIT_TO_G_PER_L[unit]
        return conc_g_L * vol_L
    elif unit in ["mmol/L", "mol/L"]:
        vol_L = (total_mass_g / density) / 1000.0
        factor = 1e-3 if unit == "mmol/L" else 1.0
        return (conc * factor) * vol_L * molar_mass
    elif unit == "% (w/w)":
        return (conc / 100.0) * total_mass_g
    elif unit == "% (v/v)":
        vol_total = total_mass_g / density
        vol_solute = (conc / 100.0) * vol_total
        return vol_solute * density 
    return 0.0

def convert_solute_to_target_unit(solute_g, total_mass_g, total_vol_ml, target_unit, ref_molar_mass):
    """换算回目标浓度单位"""
    if total_mass_g <= 1e-9 or total_vol_ml <= 1e-9: return 0.0
    
    if target_unit == "% (w/w)": return (solute_g / total_mass_g) * 100.0 
    if target_unit == "% (v/v)": return (solute_g / total_mass_g) * 100.0 

    total_vol_L = total_vol_ml / 1000.0
    if target_unit == "g/L": return solute_g / total_vol_L
    elif target_unit == "mg/L": return (solute_g * 1000) / total_vol_L
    elif target_unit == "μg/L": return (solute_g * 1e6) / total_vol_L
    elif target_unit == "mol/L": return (solute_g / ref_molar_mass) / total_vol_L
    elif target_unit == "mmol/L": return ((solute_g / ref_molar_mass) * 1000) / total_vol_L
    return 0.0

def solve_two_component_mixture(c1, d1, c2, d2, target_vol_ml, target_conc, unit):
    """解二元混合方程"""
    val_c1, val_c2, val_ct = c1, c2, target_conc
    if val_c1 == 0 and val_c2 == 0: return None, "请输入组分浓度"
    epsilon = 1e-7
    min_c, max_c = min(c1, c2), max(c1, c2)
    
    if not (min_c - epsilon <= target_conc <= max_c + epsilon):
        return None, f"目标浓度必须介于 {min_c} - {max_c} 之间"
    if abs(c1 - c2) < epsilon:
        return None, "两组分浓度相同"

    is_vol_based = unit not in ["% (w/w)"]
    if is_vol_based:
        v1 = target_vol_ml * (val_ct - val_c2) / (val_c1 - val_c2)
        v2 = target_vol_ml - v1
        return (v1 * d1, v2 * d2), None
    else:
        if abs(c1 - target_conc) < epsilon: return None, "目标浓度与组分1相同"
        ratio_m1_m2 = (target_conc - c2) / (c1 - target_conc)
        m2 = target_vol_ml / (ratio_m1_m2/d1 + 1/d2)
        m1 = m2 * ratio_m1_m2
        return (m1, m2), None

# =========================
# 4. 侧边栏：输入区域
# =========================
with st.sidebar:
    st.title("🧪 实验参数配置")
    exp_name = st.text_input("实验内容名称", value="未命名实验")

    with st.expander("🌍 环境与单位设置", expanded=True):
        room_temp = st.slider("室温 (℃)", 10.0, 35.0, 22.0, 0.5)
        d_water = get_water_density(room_temp)
        d_saline = get_saline_density(room_temp)
        st.markdown(f"<div class='ref-text'>💧 纯水密度: <b>{d_water:.4f}</b> g/mL<br>🧂 盐水密度: <b>{d_saline:.4f}</b> g/mL</div>", unsafe_allow_html=True)
        
        c_unit1, c_unit2, c_unit3 = st.columns([1.2, 1, 1])
        with c_unit1:
            # 浓度默认 mg/L (index 1)
            conc_unit = st.selectbox("浓度单位", ["μg/L", "mg/L", "g/L", "mmol/L", "mol/L", "% (w/w)", "% (v/v)"], index=1)
        with c_unit2:
            # 质量默认 mg (index 1)
            mass_unit = st.selectbox("质量单位", ["μg", "mg", "g", "kg"], index=1)
        with c_unit3:
            # 体积默认 μL (index 0)
            vol_unit = st.selectbox("体积单位", ["μL", "mL", "L"], index=0)
            
        material_count = st.number_input("混合组分数量", 2, 10, 2)

    st.markdown("---")
    
    # === 目标设置区 ===
    st.markdown("#### 🎯 目标设定 (智能计算)")
    c_tgt1, c_tgt2 = st.columns(2)
    with c_tgt1:
        target_vol_input = st.number_input(f"目标总体积 ({vol_unit})", min_value=0.0, value=0.0, step=1.0)
    with c_tgt2:
        target_conc_input = st.number_input(f"目标浓度 ({conc_unit})", min_value=0.0, value=0.0, step=0.1)
    
    target_vol_ml = target_vol_input * VOL_UNIT_TO_ML[vol_unit]
    is_auto_solve_mode = (target_vol_ml > 0 and target_conc_input > 0)
    
    if is_auto_solve_mode:
        if int(material_count) != 2:
            st.error("⚠️ 智能反算仅支持 2 种组分")
            is_valid_solve = False
        else:
            st.success(f"⚡ 实时反算模式已激活")
            is_valid_solve = True
    else:
        is_valid_solve = False

    st.markdown("---")
    st.markdown("#### 📝 组分参数录入")
    materials_input = []
    ref_molar_mass = 58.44 

    for i in range(int(material_count)):
        st.markdown(f"**🟢 组分 {i + 1}**")
        c1, c2, c3 = st.columns([1, 1, 1])
        
        # 默认值逻辑
        default_conc = 0.0
        if i == 0: default_conc = 0.0
        elif i == 1: default_conc = 100.0
        
        with c1:
            conc = st.number_input(f"浓度", min_value=0.0, value=default_conc, key=f"c_{i}")
        with c2:
            dens = st.number_input(f"密度 (g/mL)", min_value=0.1, value=1.0, step=0.001, format="%.4f", key=f"d_{i}")
        
        # 标签修改为“加入质量”
        mass_label = f"加入质量 ({mass_unit})"
        is_disabled = False
        if is_valid_solve:
            mass_label = "自动计算 (加入质量)"
            is_disabled = True
            
        with c3:
            mass = st.number_input(mass_label, min_value=0.0, value=100.0, key=f"m_{i}", disabled=is_disabled)
            
        mm = 58.44
        if conc_unit in ["mmol/L", "mol/L"]:
            mm = st.number_input(f"摩尔质量 (g/mol) - 组分{i+1}", value=58.44, key=f"mm_{i}")
            if i == 0: ref_molar_mass = mm
        
        st.markdown("<hr style='margin: 5px 0; border-top: 1px dashed #ddd;'>", unsafe_allow_html=True)
        materials_input.append({"id": i, "conc": conc, "mass": mass, "density": dens, "molar_mass": mm})

# =========================
# 5. 主逻辑计算
# =========================
final_materials = []
solve_error_msg = None

if is_valid_solve:
    m1_conc, m1_dens = materials_input[0]["conc"], materials_input[0]["density"]
    m2_conc, m2_dens = materials_input[1]["conc"], materials_input[1]["density"]
    solved_masses_g, err = solve_two_component_mixture(m1_conc, m1_dens, m2_conc, m2_dens, target_vol_ml, target_conc_input, conc_unit)
    if err:
        solve_error_msg = err
    else:
        for idx, item in enumerate(materials_input):
            calc_mass_g = solved_masses_g[idx]
            calc_vol_mL = calc_mass_g / item["density"]
            solute_g = calculate_solute_mass(item["conc"], conc_unit, item["molar_mass"], item["density"], calc_mass_g)
            final_materials.append({**item, "质量(g)": calc_mass_g, "体积(mL)": calc_vol_mL, "溶质质量(g)": solute_g})

elif target_vol_ml > 0 and not is_valid_solve:
    base_vol = sum([ (item["mass"] * MASS_UNIT_TO_G[mass_unit]) / item["density"] for item in materials_input ])
    scaling_factor = target_vol_ml / base_vol if base_vol > 0 else 0
    for item in materials_input:
        req_m_g = (item["mass"] * MASS_UNIT_TO_G[mass_unit]) * scaling_factor
        solute_g = calculate_solute_mass(item["conc"], conc_unit, item["molar_mass"], item["density"], req_m_g)
        final_materials.append({**item, "质量(g)": req_m_g, "体积(mL)": req_m_g/item["density"], "溶质质量(g)": solute_g})
else:
    for item in materials_input:
        m_g = item["mass"] * MASS_UNIT_TO_G[mass_unit]
        solute_g = calculate_solute_mass(item["conc"], conc_unit, item["molar_mass"], item["density"], m_g)
        final_materials.append({**item, "质量(g)": m_g, "体积(mL)": m_g/item["density"], "溶质质量(g)": solute_g})

df = pd.DataFrame(final_materials)

if not df.empty:
    theo_mass_g = df["质量(g)"].sum()
    theo_solute_g = df["溶质质量(g)"].sum()
    theo_vol_ml = df["体积(mL)"].sum()
    
    display_mass = theo_mass_g / MASS_UNIT_TO_G[mass_unit]
    display_vol = theo_vol_ml / VOL_UNIT_TO_ML[vol_unit]
    display_density = theo_mass_g / theo_vol_ml if theo_vol_ml > 0 else 0
    display_conc = convert_solute_to_target_unit(theo_solute_g, theo_mass_g, theo_vol_ml, conc_unit, ref_molar_mass)
else:
    display_mass, display_vol, display_density, display_conc = 0,0,0,0

# =========================
# 6. 主界面显示
# =========================
st.title(f"液体混合浓度与密度计算器 {APP_VERSION}")
st.caption(f"当前实验项目：{exp_name}")

if solve_error_msg:
    st.error(f"❌ 计算受阻: {solve_error_msg} (请在侧边栏调整参数)")
else:
    if is_valid_solve:
        st.info(f"💡 智能模式: 已根据目标浓度自动计算出各组分质量。")
    elif target_vol_ml > 0:
        st.info(f"💡 缩放模式: 已将配方缩放至目标体积。")

    # 1. 核心指标卡片
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.markdown(f'<div class="metric-card"><div class="metric-label">混合后总质量 ({mass_unit})</div><div class="metric-value">{display_mass:.2f}</div></div>', unsafe_allow_html=True)
    with c2:
        st.markdown(f'<div class="metric-card"><div class="metric-label">混合后总体积 ({vol_unit})</div><div class="metric-value">{display_vol:.2f}</div></div>', unsafe_allow_html=True)
    with c3:
        st.markdown(f'<div class="metric-card"><div class="metric-label">混合密度 (g/mL)</div><div class="metric-value">{display_density:.4f}</div></div>', unsafe_allow_html=True)
    with c4:
        st.markdown(f'<div class="metric-card"><div class="metric-label">最终浓度 ({conc_unit})</div><div class="metric-value" style="color:#d63031">{display_conc:.3f}</div></div>', unsafe_allow_html=True)

    st.divider()

    # 2. 详细配方表
    st.subheader("📋 详细配方表")
    
    display_data = []
    for _, row in df.iterrows():
        req_mass_val = row["质量(g)"] / MASS_UNIT_TO_G[mass_unit]
        solute_str = auto_format_solute(row["溶质质量(g)"])
        
        display_data.append({
            "组分名称": f"组分 {int(row['id'])+1}",
            "原始浓度": f"{row['conc']}",
            "密度 (g/mL)": f"{row['density']:.4f}",
            f"加入质量 ({mass_unit})": f"{req_mass_val:.2f}", # 修改为“加入质量”
            "含溶质 (智能单位)": solute_str
        })
    
    display_df = pd.DataFrame(display_data)
    
    # 使用 Styler 全居中
    styler = display_df.style.set_properties(**{'text-align': 'center'}) \
                             .set_table_styles([
                                 dict(selector='th', props=[('text-align', 'center')]),
                                 dict(selector='td', props=[('text-align', 'center')])
                             ])
    
    st.table(styler)

    # 3. PDF 导出
    st.divider()
    
    def footer_canvas(canvas, doc):
        canvas.saveState()
        canvas.setFont('Helvetica', 9)
        canvas.setFillColor(colors.lightgrey)
        w, h = A4
        canvas.drawRightString(w - 30, 20, f"Generated by {APP_VERSION}")
        canvas.restoreState()

    def generate_pdf():
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdf")
        doc = SimpleDocTemplate(tmp.name, pagesize=A4)
        frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
        template = PageTemplate(id='test', frames=frame, onPage=footer_canvas)
        doc.addPageTemplates([template])

        styles = getSampleStyleSheet()
        style_title = ParagraphStyle('TitleCN', parent=styles['Title'], fontName=FONT_NAME, fontSize=22, spaceAfter=20)
        style_h2 = ParagraphStyle('H2CN', parent=styles['Heading2'], fontName=FONT_NAME, fontSize=14, spaceBefore=15, spaceAfter=10)
        style_normal = ParagraphStyle('NormalCN', parent=styles['Normal'], fontName=FONT_NAME, fontSize=10, leading=14, alignment=1) # 1=Center
        style_left = ParagraphStyle('LeftCN', parent=styles['Normal'], fontName=FONT_NAME, fontSize=10, leading=14)

        elements = []
        elements.append(Paragraph(f"实验报告：{exp_name}", style_title))
        elements.append(Paragraph(f"生成时间: {time.strftime('%Y-%m-%d %H:%M')}", style_normal))
        elements.append(Spacer(1, 15))
        
        elements.append(Paragraph("1. 环境与目标", style_h2))
        env_text = f"""
        <b>室温:</b> {room_temp} ℃ <br/>
        <b>目标单位:</b> {conc_unit} (浓度) | {vol_unit} (体积)<br/>
        <b>参考密度:</b> 纯水 ({d_water:.4f} g/mL) | 生理盐水 ({d_saline:.4f} g/mL)
        """
        elements.append(Paragraph(env_text, style_left))
        
        elements.append(Paragraph("2. 混合结果总览", style_h2))
        res_text = f"""
        <b>总质量:</b> {display_mass:.2f} {mass_unit}<br/>
        <b>总体积:</b> {display_vol:.2f} {vol_unit}<br/>
        <b>混合密度:</b> {display_density:.4f} g/mL<br/>
        <b>混合浓度:</b> {display_conc:.3f} {conc_unit}
        """
        elements.append(Paragraph(res_text, style_left))
        
        elements.append(Paragraph("3. 详细配方表", style_h2))
        # 修改PDF表头为“加入质量”
        headers = ["组分", f"原始浓度", f"密度\n(g/mL)", f"加入质量\n({mass_unit})", "含溶质\n(自动单位)"]
        data = [headers]
        for _, row in df.iterrows():
            req_mass = row["质量(g)"] / MASS_UNIT_TO_G[mass_unit]
            solute_str = auto_format_solute(row["溶质质量(g)"])
            data.append([
                f"组分 {int(row['id'])+1}",
                f"{row['conc']}",
                f"{row['density']:.4f}",
                f"{req_mass:.2f}",
                solute_str
            ])
            
        t = Table(data, colWidths=[90, 80, 80, 90, 100])
        t.setStyle(TableStyle([
            ('FONTNAME', (0,0), (-1,-1), FONT_NAME),
            ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
            ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#e6e6e6")),
            ('ALIGN', (0,0), (-1,-1), 'CENTER'),
            ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ]))
        elements.append(t)
        doc.build(elements)
        return tmp.name

    col_btn, col_empty = st.columns([1, 4])
    with col_btn:
        if st.button("📥 生成 PDF 报告", type="primary"):
            try:
                path = generate_pdf()
                with open(path, "rb") as f:
                    st.download_button(f"下载PDF报告", data=f, file_name=f"{exp_name}.pdf", mime="application/pdf")
            except Exception as e:
                st.error(f"PDF生成错误: {e}")

