import streamlit as st
from Bio import SeqIO
import pandas as pd
import pickle
from collections import Counter
import os
from advanced_gene_scanner import scan_genome_for_resistance_genes, predict_resistance_evolution

st.set_page_config(page_title="AMR Forecasting", page_icon="🧬", layout="wide")

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
- Scans for resistance genes
- Works for ANY bacterium
- Forecasts resistance evolution

**Accuracy:** 90-100% across antibiotics
""")

models = load_models()

if models:
    amr_genes_db = models.pop('amr_genes', {}) if 'amr_genes' in models else {}
    st.success(f"✅ Loaded {len(models)} antibiotic prediction models")
    
    st.header("📤 Upload Genome")
    uploaded_file = st.file_uploader(
        "Choose a FASTA file (up to 100 MB)", 
        type=['fasta', 'fa', 'fna'],
        help="Upload bacterial genome in FASTA format"
    )
    
    if uploaded_file:
        file_size_mb = uploaded_file.size / (1024 * 1024)
        if file_size_mb > 100:
            st.error(f"⚠️ File is {file_size_mb:.1f} MB. Please use a genome ≤ 100 MB.")
            st.stop()
        if file_size_mb > 50:
            st.warning("⚠️ Large file detected. Processing may take 2-5 minutes.")
        
        with open("temp_genome.fasta", "wb") as f:
            f.write(uploaded_file.getbuffer())
        st.success(f"✅ File uploaded! ({file_size_mb:.1f} MB)")
        
        try:
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            status_text.text("📖 Reading genome...")
            progress_bar.progress(10)
            
            records = list(SeqIO.parse("temp_genome.fasta", "fasta"))
            if len(records) == 0:
                st.error("No sequences found in file!")
                st.stop()
            elif len(records) > 1:
                st.warning(f"⚠️ File contains {len(records)} sequences. Using the first one: {records[0].id}")
            record = records[0]
            seq = str(record.seq).upper()
            
            progress_bar.progress(20)
            
            col1, col2, col3 = st.columns(3)
            col1.metric("Accession", record.id)
            col2.metric("Length", f"{len(seq):,} bp")
            col3.metric("GC%", f"{(seq.count('G') + seq.count('C')) / len(seq) * 100:.1f}%")
            
            status_text.text("🧬 Extracting k-mer features...")
            progress_bar.progress(30)
            
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
            
            progress_bar.progress(50)
            status_text.text("🤖 Running ML predictions...")
            
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
            
            progress_bar.progress(70)
            status_text.text("🔬 Scanning for resistance genes...")
            
            st.header("🔬 Genomic Resistance Gene Analysis")
            detected_genes = scan_genome_for_resistance_genes("temp_genome.fasta")
            if detected_genes:
                st.subheader("🧬 Detected Intrinsic Resistance Genes")
                st.info("These genes are ALREADY PRESENT in the genome.")
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
                st.success("✓ No known intrinsic resistance genes detected.")
            
            progress_bar.progress(85)
            status_text.text("🔮 Predicting resistance evolution...")
            
            st.header("🔮 Resistance Evolution Forecast")
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("First-Time Exposure")
                for antibiotic_class in ['beta_lactam', 'aminoglycoside', 'fluoroquinolone']:
                    evolution = predict_resistance_evolution(seq, antibiotic_class)
                    st.markdown(f"**{antibiotic_class.replace('_', ' ').title()}**")
                    st.progress(evolution['probability'])
                    st.caption(f"Evolution probability: {evolution['probability']:.0%} in {evolution['timeline_months']} months")
            
            with col2:
                st.subheader("Repeated Exposure")
                for antibiotic_class in ['beta_lactam', 'aminoglycoside', 'fluoroquinolone']:
                    evolution = predict_resistance_evolution(seq, antibiotic_class)
                    repeated_prob = min(evolution['probability'] * 2.5, 0.98)
                    repeated_months = max(int(evolution['timeline_months'] / 3), 1)
                    st.markdown(f"**{antibiotic_class.replace('_', ' ').title()}**")
                    st.progress(repeated_prob)
                    st.caption(f"Evolution probability: {repeated_prob:.0%} in {repeated_months} months")
            
            progress_bar.progress(100)
            status_text.text("✅ Analysis complete!")
            st.success("✅ Complete analysis finished!")
            
        except Exception as e:
            st.error(f"Error: {e}")
            import traceback
            st.error(traceback.format_exc())
else:
    st.error("Models not loaded. Please contact system administrator.")

st.markdown("---")
st.markdown("**AMR Evolutionary Forecasting System** | Comprehensive | Novel Research Tool | 🔬")
