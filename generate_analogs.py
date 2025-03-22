import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, Crippen, rdMolDescriptors, QED
import pandas as pd
import random
import base64
from io import BytesIO

# ---------------------- Configuration ----------------------
st.set_page_config(page_title="Compound Analog Generator", page_icon="🧪", layout="wide")
st.title("🧪 Compound Analog Generator")
st.markdown("Generate modified analogs with molecular properties from SMILES input")

# ---------------------- Core Functions ----------------------
def calculate_properties(smiles):
    """Calculate molecular properties and drug-likeness metrics"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        logp = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        return {
            'MW': round(Descriptors.MolWt(mol), 2),
            'LogP': round(logp, 2),
            'TPSA': round(tpsa, 2),
            'HBD': rdMolDescriptors.CalcNumHBD(mol),
            'HBA': rdMolDescriptors.CalcNumHBA(mol),
            'Rot Bonds': Descriptors.NumRotatableBonds(mol),
            'Lipinski': "Yes" if sum([
                Descriptors.MolWt(mol) > 500, 
                logp > 5,
                rdMolDescriptors.CalcNumHBD(mol) > 5
            ]) <= 1 else "No",
            'QED': round(QED.qed(mol), 2),
            'BBB Penetration': "Yes" if (1 <= logp <= 4 and tpsa < 90) else "No"
        }
    return None

def generate_analogs(base_smiles, num_analogs=10):
    """Generate molecular analogs"""
    analogs = []
    base_mol = Chem.MolFromSmiles(base_smiles)
    
    if base_mol:
        for _ in range(num_analogs):
            try:
                new_mol = Chem.RWMol(base_mol)
                atom_idx = random.randint(0, new_mol.GetNumAtoms()-1)
                new_atom = random.choice([6, 7, 8])  # C, N, O
                new_mol.AddAtom(Chem.Atom(new_atom))
                new_mol.AddBond(atom_idx, new_mol.GetNumAtoms()-1, Chem.BondType.SINGLE)
                
                analog_smiles = Chem.MolToSmiles(new_mol)
                analogs.append({
                    'smiles': analog_smiles,
                    'props': calculate_properties(analog_smiles)
                })
            except Exception as e:
                st.error(f"Error generating analog: {str(e)}")
    return analogs

# ---------------------- UI Components ----------------------
def main():
    # Parent Compound Section
    st.sidebar.header("⚗️ Parent Compound")
    smiles_input = st.sidebar.text_input("Enter SMILES:", value="CCO")
    
    parent_props = None
    parent_mol = Chem.MolFromSmiles(smiles_input)
    if parent_mol:
        buf = BytesIO()
        Draw.MolToImage(parent_mol).save(buf, format='PNG')
        st.sidebar.image(buf.getvalue(), caption="Parent Structure", use_container_width=True)
        parent_props = calculate_properties(smiles_input)
        if parent_props:
            st.sidebar.markdown("### Parent Properties")
            parent_df = pd.DataFrame([parent_props]).T.rename(columns={0: 'Value'})
            st.sidebar.table(parent_df)
    else:
        st.sidebar.error("Invalid SMILES input")

    # Analog Generation Controls
    st.sidebar.header("⚙️ Settings")
    num_analogs = st.sidebar.slider("Number of analogs", 5, 50, 15)
    qed_threshold = st.sidebar.slider("QED Threshold", 0.0, 1.0, 0.18)
    
    if st.sidebar.button("Generate Analogs"):
        st.session_state.analogs = None
        
        if parent_mol and parent_props:
            analogs = generate_analogs(smiles_input, num_analogs)
            valid_analogs = [a for a in analogs if a['props']]
            
            if valid_analogs:
                # Create DataFrame with images
                df = pd.DataFrame([{
                    **a['props'], 
                    'SMILES': a['smiles'],
                    'Structure': Draw.MolToImage(Chem.MolFromSmiles(a['smiles']))
                } for a in valid_analogs])
                
                df.insert(0, 'Sr.No', range(1, len(df)+1))
                st.session_state.analogs = df
                st.session_state.qed_threshold = qed_threshold

    # Display results
    if 'analogs' in st.session_state and st.session_state.analogs is not None:
        filtered_df = st.session_state.analogs[
            st.session_state.analogs['QED'] >= st.session_state.qed_threshold
        ].reset_index(drop=True)
        filtered_df['Sr.No'] = range(1, len(filtered_df)+1)

        st.header(f"🔬 Generated Analogs ({len(filtered_df)} shown/{num_analogs} generated)")
        
        # Property Distribution Charts
        st.markdown("### Property Analysis")
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**Molecular Weight vs LogP**")
            st.line_chart(filtered_df.set_index('MW')['LogP'])
        with col2:
            st.markdown("**Drug-likeness Distribution (QED)**")
            st.bar_chart(filtered_df['QED'])
        
        # Molecular Structures Grid
        st.markdown("### Analog Structures")
        cols = st.columns(4)
        for idx, row in filtered_df.iterrows():
            with cols[idx % 4]:
                st.image(row['Structure'], 
                       caption=f"Sr.No: {row['Sr.No']}\nQED: {row['QED']:.2f} | BBB: {row['BBB Penetration']}",
                       use_container_width=True)
        
        # Data Table
        st.markdown("### Detailed Properties")
        st.dataframe(filtered_df.drop(columns=['Structure']).set_index('Sr.No'))

        # CSV Download
        csv = filtered_df.drop(columns=['Structure']).to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()
        st.markdown(f'<a href="data:file/csv;base64,{b64}" download="analogs.csv">📥 Download Analog Data</a>', 
                    unsafe_allow_html=True)

if __name__ == "__main__":
    main()