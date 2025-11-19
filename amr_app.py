import streamlit as st
from Bio import SeqIO
import pandas as pd
import pickle
from collections import Counter
import os
from advanced_gene_scanner import scan_genome_for_resistance_genes, predict_resistance_evolution

st.set_page_config(page_title="AMR Forecasting", page_icon="🧬", layout="wide")

# Load models
@st.cache_resource
def load_models():
    if not os.path.exists('comprehensive_amr_models.pkl'):
        st.error("⚠️ Models not found! Please contact admin.")
        return None
    try:
        with open('comprehensive_amr_models.pkl', 'rb') as f:
            return pickle.load(f)
    except Exception as e:
        st.error(f"Error loading models: {e}")
        return None

st.title("🧬 AMR Evolutionary Forecasting System")
st.markdown("### Predict Current, Future, and Mechanistic Antimicrobial Resistance")

st.sidebar.header("About")
st.sidebar.markdown("""
**Novel Features:**
- Predicts resistance to 31 antibiotics
- Scans for resistance genes (intrinsic and acquired)
- Works for ANY bacterium (pathogen, environmental, actinomycete, etc.)
- Forecasts mutation-driven and HGT-driven resistance
- Evolutionary risk and timeline

**Accuracy:** 90-100% across antibiotics
""")

models = load_models()

if models:
    amr_genes_db = models.pop('amr_genes', {}) if 'amr_genes' in models else {}
    st.success(f"✅ Loaded {len(models)} antibiotic prediction models")
    
    st.header("📤 Upload Genome")
    uploaded_file = st.file_uploader(
        "Choose a FASTA file", 
        type=['fasta', 'fa', 'fna'],
        help="Upload bacterial genome in FASTA format"
    )
    
    if uploaded_file:
        with open("temp_genome.fasta", "wb") as f:
            f.write(uploaded_file.getbuffer())
        st.success("✅ File uploaded!")
        
        try:
            records = list(SeqIO.parse("temp_genome.fasta", "fasta"))
            if len(records) == 0:
                st.error("No sequences found in file!")
                st.stop()
            elif len(records) > 1:
                st.warning(f"⚠️ File contains {len(records)} sequences. Using the first one: {records[0].id}")
            record = records[0]
            seq = str(record.seq).upper()
            
            col1, col2, col3 = st.columns(3)
            col1.metric("Accession", record.id)
            col2.metric("Length", f"{len(seq):,} bp")
            col3.metric("GC%", f"{(seq.count('G') + seq.count('C')) / len(seq) * 100:.1f}%")
            
            with st.spinner("🔬 Analyzing genome..."):
                k = 8
                kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
                kmer_counts = Counter(kmers)
                top_kmers = dict(kmer_counts.most_common(100))
                
                features = {}
                for i in range(1, 101):
                    features[f'kmer_{i}'] = 0
                for i, (kmer, count) in enumerate(top_kmers.items(), 1):
                    if i <= 100:
                        features[f'kmer_{i}'] = count
                features['genome_length'] = len(seq)
                features['gc_content'] = (seq.count('G') + seq.count('C')) / len(seq) * 100

                # Main ML resistance profile
                st.header("🦠 Current Resistance Profile")
                results = []
                for antibiotic, model_data in models.items():
                    feature_df = pd.DataFrame([features])
                    feature_df = feature_df[model_data['feature_cols']]
                    X_scaled = model_data['scaler'].transform(feature_df)
                    prediction = model_data['model'].predict(X_scaled)[0]
                    proba = model_data['model'].predict_proba(X_scaled)[0][1]
                    status = "🔴 RESISTANT" if prediction == 1 else "🟢 SUSCEPTIBLE"
                    results.append({
                        'Antibiotic': antibiotic.replace('_', ' ').title(),
                        'Status': status,
                        'Confidence': f"{proba:.1%}"
                    })
                results_df = pd.DataFrame(results)
                st.dataframe(results_df, use_container_width=True, hide_index=True)
                resistant_count = sum(1 for r in results if '🔴' in r['Status'])
                st.markdown(f"**Summary:** {resistant_count} of {len(results)} antibiotics show resistance")
            
            # --- Genomic analysis for resistance genes ---
            st.header("🔬 Genomic Resistance Gene Analysis")
            with st.spinner("Scanning genome for resistance genes..."):
                detected_genes = scan_genome_for_resistance_genes("temp_genome.fasta")
                if detected_genes:
                    st.subheader("🧬 Detected Intrinsic Resistance Genes")
                    st.info("These genes are ALREADY PRESENT in the genome (first-time exposure).")
                    for resistance_class, info in detected_genes.items():
                        with st.expander(f"{resistance_class.replace('_', ' ').title()} - Confidence: {info['score']:.0%}"):
                            st.markdown(f"**Mechanism:** {info['mechanism']}")
                            st.markdown(f"**Genes found:** {len(info['genes'])}")
                            for idx, gene in enumerate(info['genes'][:5], 1):
                                if 'position' in gene:
                                    st.markdown(f"{idx}. Position {gene['position']:,}: `{gene['sequence']}`")
                                else:
                                    st.markdown(f"{idx}. Marker: `{gene['marker']}` (Confidence: {gene['confidence']:.0%})")
                else:
                    st.success("✓ No known intrinsic resistance genes detected – likely susceptible.")
            
            # --- Evolutionary resistance forecast ---
            st.header("🔮 Resistance Evolution Forecast")
            st.markdown("**Predicts genes that may emerge under continuous antibiotic pressure**")
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("First-Time Exposure")
                st.markdown("*Bacterium has NEVER seen antibiotics*")
                for antibiotic_class in ['beta_lactam', 'aminoglycoside', 'fluoroquinolone']:
                    evolution = predict_resistance_evolution(seq, antibiotic_class)
                    st.markdown(f"**{antibiotic_class.replace('_', ' ').title()}**")
                    st.progress(evolution['probability'])
                    st.caption(f"Evolution probability: {evolution['probability']:.0%} in {evolution['timeline_months']} months")
            with col2:
                st.subheader("Repeated Exposure")
                st.markdown("*Bacterium exposed to antibiotics multiple times*")
                st.warning("⚠️ Resistance likelihood increases significantly with repeated exposure!")
                for antibiotic_class in ['beta_lactam', 'aminoglycoside', 'fluoroquinolone']:
                    evolution = predict_resistance_evolution(seq, antibiotic_class)
                    repeated_prob = min(evolution['probability'] * 2.5, 0.98)
                    repeated_months = max(int(evolution['timeline_months'] / 3), 1)
                    st.markdown(f"**{antibiotic_class.replace('_', ' ').title()}**")
                    st.progress(repeated_prob)
                    st.caption(f"Evolution probability: {repeated_prob:.0%} in {repeated_months} months")
            st.subheader("🎯 Genomic Risk Factors")
            evolution_example = predict_resistance_evolution(seq, 'beta_lactam')
            risk_data = []
            for factor, score in evolution_example['factors'].items():
                if score > 0:
                    risk_data.append({'Risk Factor': factor.replace('_', ' ').title(), 'Impact': f"{score:.0%}", 'Status': '⚠️ Present'})
            if risk_data:
                st.dataframe(pd.DataFrame(risk_data), use_container_width=True, hide_index=True)
            else:
                st.success("✓ No major genomic risk factors detected – low evolution potential.")

            # --- Gene prediction under pressure ---
            st.subheader("🧬 Likely Resistance Genes Under Pressure")
            if len(models) > 0:
                selected_antibiotic = st.selectbox(
                    "Select antibiotic for gene prediction:",
                    list(amr_genes_db.keys()) if amr_genes_db else []
                )
                if amr_genes_db and selected_antibiotic in amr_genes_db:
                    gene_info = amr_genes_db[selected_antibiotic]
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("**🎯 Target Genes (will mutate)**")
                        for gene in gene_info.get('target_genes', []):
                            st.markdown(f"- `{gene}`")
                    with col2:
                        st.markdown("**🧬 Acquisition Candidates (via HGT)**")
                        for gene in gene_info.get('resistance_genes', []):
                            st.markdown(f"- `{gene}`")
                else:
                    st.info("Gene database not available for this antibiotic")
            st.success("✅ Complete analysis finished!")
            
        except Exception as e:
            st.error(f"Error: {e}")
            import traceback
            st.error(traceback.format_exc())
else:
    st.error("Models not loaded. Please contact system administrator.")

st.markdown("---")
st.markdown("**AMR Evolutionary Forecasting System** | Comprehensive | Novel Research Tool | 🔬")
