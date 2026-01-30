import streamlit as st
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application

# --- 1. CONFIG & STYLING ---
st.set_page_config(page_title="2BAC Math Solver - SBAKSI", page_icon="📐", layout="wide")

st.markdown("""
    <style>
    .main { background-color: #0f172a; color: #f1f5f9; }
    .math-card {
        background: rgba(30, 41, 59, 0.7);
        padding: 24px;
        border-radius: 15px;
        border: 1px solid rgba(56, 189, 248, 0.3);
        margin-bottom: 20px;
    }
    .redaction-box {
        background: #1e293b;
        padding: 15px;
        border-radius: 10px;
        border-left: 5px solid #10b981;
        color: #e2e8f0;
        margin: 10px 0;
        font-size: 1rem;
    }
    .footer {
        position: fixed;
        left: 0;
        bottom: 0;
        width: 100%;
        background-color: #1e293b;
        color: #38bdf8;
        text-align: center;
        padding: 8px;
        font-weight: bold;
        border-top: 1px solid #38bdf8;
        z-index: 100;
    }
    </style>
    """, unsafe_allow_html=True)

# --- 2. ADVANCED MATH ENGINE ---

def format_moroccan_interval(domain_set):
    if domain_set.is_EmptySet: return "∅"
    if domain_set == sp.S.Reals: return "ℝ"
    def process(interval):
        if isinstance(interval, sp.FiniteSet):
            return "ℝ \\ {" + ", ".join([sp.latex(x) for x in interval.args]) + "}"
        if isinstance(interval, sp.Interval):
            left = "]" if interval.left_open else "["
            right = "[" if interval.right_open else "]"
            start = sp.latex(interval.start).replace(r"\infty", "+∞").replace(r"-\infty", "-∞")
            end = sp.latex(interval.end).replace(r"\infty", "+∞").replace(r"-\infty", "-∞")
            return f"{left}{start} ; {end}{right}"
        return str(interval)
    return " ∪ ".join([process(i) for i in (domain_set.args if isinstance(domain_set, sp.Union) else [domain_set])])

def get_accurate_merged_redaction(f, f_prime, x, domain):
    """The Final Fixed Logic for Merged Closed/Open Intervals."""
    sentences = []
    try:
        # Get every point where the sign could change or function could stop
        zeros = sp.solve(f_prime, x)
        singularities = sp.singularities(f, x)
        crit_pts = sorted([float(p.evalf()) for p in (list(zeros) + list(singularities)) if p.is_real])
        
        bounds = [-float('inf')] + crit_pts + [float('inf')]
        raw_segments = []
        
        for i in range(len(bounds)-1):
            a, b = bounds[i], bounds[i+1]
            mid = (a + b) / 2 if abs(a) < 1e6 and abs(b) < 1e6 else (a+1 if b > 1e6 else b-1)
            
            if domain.contains(mid):
                sign_val = f_prime.subs(x, mid).evalf()
                sign = "pos" if sign_val > 0 else "neg"
                raw_segments.append({'start': a, 'end': b, 'sign': sign})

        if not raw_segments: return ["Aucune variation détectée."]

        # Merge segments with the same sign
        merged = []
        curr = raw_segments[0]
        for next_seg in raw_segments[1:]:
            if next_seg['sign'] == curr['sign']:
                curr['end'] = next_seg['end']
            else:
                merged.append(curr)
                curr = next_seg
        merged.append(curr)

        # Apply Moroccan Bracket Rules (Check continuity at boundaries)
        for m in merged:
            l_bracket, r_bracket = "]", "["
            l_val = "-∞" if m['start'] == -float('inf') else round(m['start'], 2)
            r_val = "+∞" if m['end'] == float('inf') else round(m['end'], 2)
            
            # Close bracket if point belongs to Domain
            if m['start'] != -float('inf') and domain.contains(m['start']): l_bracket = "["
            if m['end'] != float('inf') and domain.contains(m['end']): r_bracket = "]"
            
            var_text = "**strictement croissante** ↗️" if m['sign'] == "pos" else "**strictement décroissante** ↘️"
            sentences.append(f"Sur {l_bracket}{l_val} ; {r_val}{r_bracket}, $f$ est {var_text}")
            
    except: return ["Calcul des intervalles en cours..."]
    return sentences

def get_geometric_interpretation(f, x, df):
    items = []
    # Vertical Asymptotes
    sings = sp.singularities(f, x)
    for s in sings:
        try:
            if s.is_real:
                if sp.limit(f, x, s) in [sp.oo, -sp.oo]:
                    items.append(f"✅ **Asymptote Verticale :** $x = {sp.latex(s)}$")
        except: pass
    # Branches at Infinity
    for d, label in [(sp.oo, "+∞"), (-sp.oo, "-∞")]:
        try:
            L = sp.limit(f, x, d)
            if L.is_finite: items.append(f"✅ **Asymptote Horizontale :** $y = {sp.latex(L)}$ au voisinage de {label}")
            elif L in [sp.oo, -sp.oo]:
                a = sp.limit(f/x, x, d)
                if a == 0: items.append(f"✅ **Branche Parabolique :** Direction $(Ox)$ au voisinage de {label}")
                elif a in [sp.oo, -sp.oo]: items.append(f"✅ **Branche Parabolique :** Direction $(Oy)$ au voisinage de {label}")
                else:
                    b = sp.limit(f - a*x, x, d)
                    if b.is_finite: items.append(f"✅ **Asymptote Oblique :** $y = {sp.latex(a)}x + {sp.latex(b)}$ au voisinage de {label}")
                    else: items.append(f"✅ **Branche Parabolique :** Direction $y = {sp.latex(a)}x$ au voisinage de {label}")
        except: pass
    return items

# --- 3. UI LAYOUT ---
with st.sidebar:
    st.markdown("## ⚙️ SBAKSI Menu")
    if st.button("📐 Racine: sqrt(x)"): st.session_state.f_in = "sqrt(x)"
    if st.button("📈 ln(x)/x"): st.session_state.f_in = "ln(x)/x"
    if st.button("📉 (x+1)/(x-1)"): st.session_state.f_in = "((x+1)/(x-1))"
    st.divider()
    st.markdown("### Rappel Racine\n`sqrt(x)` ou `x**(1/n)`")

st.title("🎓 SBAKSI Math 2BAC Suite")

if 'f_in' not in st.session_state: st.session_state.f_in = "sqrt(x)"
user_input = st.text_input("Fonction f(x) :", value=st.session_state.f_in)

if user_input:
    try:
        clean = user_input.replace('racine', 'sqrt').replace('√', 'sqrt').replace('ln', 'log')
        trans = standard_transformations + (implicit_multiplication_application,)
        x = sp.symbols('x')
        f = parse_expr(clean, transformations=trans)
        f_prime = sp.simplify(sp.diff(f, x))
        df = sp.calculus.util.continuous_domain(f, x, sp.S.Reals)

        st.markdown(f'<div class="math-card" style="text-align:center;"><h1>$f(x) = {sp.latex(f)}$</h1></div>', unsafe_allow_html=True)

        t1, t2, t3 = st.tabs(["📌 Etude & Rédaction", "🔍 Interprétation Géo", "📊 Construction"])

        with t1:
            c1, c2 = st.columns([1, 1.2])
            with c1:
                st.subheader("1. Domaine et Dérivée")
                st.info(f"**$D_f = {format_moroccan_interval(df)}$**")
                st.latex(rf"f'(x) = {sp.latex(f_prime)}")
                
                

            with c2:
                st.subheader("2. Rédaction (Intervalles Fixés)")
                lines = get_accurate_merged_redaction(f, f_prime, x, df)
                for line in lines:
                    st.markdown(f"<div class='redaction-box'>{line}</div>", unsafe_allow_html=True)

        with t2:
            st.subheader("🚩 Analyse des Branches & Asymptotes")
            geo = get_geometric_interpretation(f, x, df)
            for g in geo: st.markdown(f"<div class='redaction-box'>{g}</div>", unsafe_allow_html=True)
            
            

        with t3:
            st.subheader("🛠️ Construction (Focus Racine)")
            st.markdown("- **Point d'arrêt :** Si $D_f = [a, ...$, placez le point $(a, f(a))$ en premier.\n- **Tangente :** Si $f'(a)$ est infinie, dessinez une demi-tangente verticale.")
            f_n = sp.lambdify(x, f, modules=['numpy'])
            x_v = np.linspace(-10, 10, 1000)
            with np.errstate(all='ignore'):
                y_v = f_n(x_v)
                y_v[np.abs(y_v) > 30] = np.nan
            fig, ax = plt.subplots(facecolor='#1e293b')
            ax.set_facecolor('#0f172a')
            ax.plot(x_v, y_v, color='#38bdf8', lw=2)
            ax.axhline(0, color='white', alpha=0.3); ax.axvline(0, color='white', alpha=0.3)
            st.pyplot(fig)

    except Exception as e:
        st.error(f"Saisie invalide : {e}")

st.markdown(f'<div class="footer">Made by SBAKSI 🇲🇦</div>', unsafe_allow_html=True)