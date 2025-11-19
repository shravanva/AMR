import streamlit as st
import os

st.title("🧬 AMR Test")

# Check if model file exists
if os.path.exists('comprehensive_amr_models.pkl'):
    file_size = os.path.getsize('comprehensive_amr_models.pkl') / (1024 * 1024)
    st.success(f"✅ Model file found! Size: {file_size:.1f} MB")
else:
    st.error("❌ Model file NOT found!")

st.write("If you see this, the app is working!")
