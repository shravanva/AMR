import streamlit as st
import os

st.title("AMR App Diagnostic")

if os.path.exists('comprehensive_amr_models.pkl'):
    st.success("Found comprehensive_amr_models.pkl!")
    st.write(f"Size: {os.path.getsize('comprehensive_amr_models.pkl') / (1024*1024):.1f} MB")
else:
    st.error("comprehensive_amr_models.pkl NOT FOUND in deployment.")

st.write("If you see this message (even without the model), the app runs. If it just hangs, Cloud can't see your model file!")
