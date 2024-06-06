import streamlit as st
import joblib
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import base64

# Title of the app with customized heading style
st.markdown("<h1 style='color: black; font-size: 36px; font-weight: bold;'>Protein Sequence Alignment</h1>", unsafe_allow_html=True)

# Load the trained model using joblib
model = joblib.load(r'C:\Users\megha\Documents\Streamlit\best_model.joblib')

# Function to encode image to base64
def encode_image(image_file):
    with open(image_file, "rb") as f:
        encoded_string = base64.b64encode(f.read()).decode()
    return encoded_string

# Custom CSS for styling including background image
background_image = encode_image('protein image.jpg')  # Replace with the filename of your background image
st.markdown(f"""
    <style>
        .stApp {{
            background-image: url('data:image/jpg;base64,{background_image}');
            background-size: cover;
            background-repeat: no-repeat;
            color: black;
        }}
        .stTextInput {{
            font-size: 14px;
            font-weight: bold;
            color: #1E90FF;
            height: 150px !important;
        }}
        .stButton button {{
            background-color: #4CAF50;
            color: white;
        }}
        .title {{
            font-size: 36px;
            font-weight: bold;
            color: black;
        }}
        .result {{
            font-size: 20px;
            font-weight: bold;
            color: #FF5722;
        }}
    </style>
""", unsafe_allow_html=True)

# Function to analyze sequence features
def analyze_sequence(sequence):
    analysis = ProteinAnalysis(sequence)
    return [
        analysis.gravy(),
        analysis.isoelectric_point(),
        analysis.aromaticity(),
        analysis.instability_index(),
        analysis.molecular_weight(),
        *analysis.secondary_structure_fraction()
    ]

# Function to clean sequence
def clean_sequence(sequence):
    valid_amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join([aa for aa in sequence.upper() if aa in valid_amino_acids])

# Input fields for protein sequences with customized styling
st.markdown("**Protein A Sequence**")
protein_a_input = st.text_area("", 
                               placeholder="Enter Protein A sequence here...", 
                               help="Paste your Protein A sequence.",
                               height=100, 
                               max_chars=10000,
                               key='protein_a_input')
st.markdown("**Protein B Sequence**")
protein_b_input = st.text_area("", 
                               placeholder="Enter Protein B sequence here...", 
                               help="Paste your Protein B sequence.",
                               height=100, 
                               max_chars=10000,
                               key='protein_b_input')

# Button to predict identity score
if st.button("Predict Identity Score"):
    with st.spinner("Analyzing sequences and predicting..."):
        if protein_a_input and protein_b_input:
            # Clean and analyze sequences
            protein_a_cleaned = clean_sequence(protein_a_input)
            protein_b_cleaned = clean_sequence(protein_b_input)
            
            features_a = analyze_sequence(protein_a_cleaned)
            features_b = analyze_sequence(protein_b_cleaned)
            
            # Create a DataFrame for the model
            feature_columns = ['gravy', 'isoelectric_point', 'aromaticity', 'instability_index', 'molecular_weight', 'helix', 'turn', 'sheet']
            data = {}
            
            for i, col in enumerate(feature_columns):
                data[f'{col}_A'] = [features_a[i]]
                data[f'{col}_B'] = [features_b[i]]
            
            df = pd.DataFrame(data)
            
            # Predict identity score
            identity_score = model.predict(df)[0]
            
            # Display the result
            st.markdown(f"<div class='result'>**Predicted Identity Score: {identity_score:.2f}%**</div>", unsafe_allow_html=True)
        else:
            st.warning("Please input both protein sequences.")
