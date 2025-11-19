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
- Scans for resistance genes
- Works for ANY bacterium
- Forecasts resistance evolution

**Accuracy:** 90-100% across antibiotics

**Note:** For optimal performance, use genomes < 10 MB
""")

models = load_models()

if models:
    amr_genes_db = models.pop('amr_genes', {}) if 'amr_genes' in models else {}
    st.success(f"✅ Loaded {len(models)} antibiotic prediction models")
    
    st.header("📤 Upload Genome")
    uploaded_file = st.file_uploader(
        "Choose a FASTA file (recommended < 10 MB)", 
        type=['fasta', 'fa', 'fna'],
        help="Upload bacterial genome in FASTA format"
    )
    
    if uploaded_file:
        # Check file size
        file_size_mb = uploaded_file.size / (1024 * 1024)
        if file_size_mb > 20:
            st.error(f"⚠️ File too large ({file_size_mb:.1f} MB). Please use a genome < 20 MB.")
            st.stop()
        
        with open("temp_genome.fasta", "wb") as f:
            f.write(uploaded_file.getbuffer())
        st.success(f"✅ File uploaded! ({file_size_mb:.1f} MB)")
        
        try:
            # Progress bar
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            status_text.text("Reading genome...")
            progress_bar.progress(10)
            
            records = list(SeqIO.parse("temp_genome.fasta", "fasta"))
            if len(records) == 0:
                st.error("No sequences found in file!")
                st.stop()
            elif len(records) > 1:
                st.warning(f"
