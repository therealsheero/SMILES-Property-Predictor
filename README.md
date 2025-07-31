
# 🧬 SMILES Property Predictor

A Python-based tool to predict **molecular properties** directly from **SMILES (Simplified Molecular Input Line Entry System)** strings using a combination of **RDKit** (rule-based cheminformatics) and **Machine Learning models** for advanced predictions like **solubility (LogS)** and **toxicity**.

---

## 🧪 What Are SMILES?

**SMILES** (Simplified Molecular Input Line Entry System) is a compact string notation used to describe the structure of chemical molecules. For example:
```

Aspirin → `CC(=O)Oc1ccccc1C(=O)O`

````

These strings encode atoms, bonds, rings, and branches in a form that can be easily parsed and analyzed by cheminformatics tools like RDKit.

---

## 🧠 What This Project Does

Given a SMILES string, this tool predicts **12 important molecular properties**, including:

|  No. | Property              |  Description                                      | Why It's Important                                                              |
| -----: | ------------------------ | ---------------------------------------------------- | ------------------------------------------------------------------------------- |
|      1 | **Canonical SMILES**     | A standardized form of the SMILES string.            | Helps ensure consistency when comparing or storing molecules.                   |
|      2 | **Molecular Formula**    | Elemental composition (e.g. C9H8O4).                 | Fundamental chemical identity — needed in molecular databases and search.       |
|      3 | **Molecular Weight**     | Total molecular mass in Daltons (g/mol).             | Useful for pharmacokinetics and dosing calculations.                            |
|      4 | **LogP**                 | Octanol–water partition coefficient (lipophilicity). | Measures how hydrophobic a compound is — affects absorption, distribution, etc. |
|      5 | **H-Bond Donors**        | Count of -OH or -NH groups.                          | Affects solubility and interaction with biological targets.                     |
|      6 | **H-Bond Acceptors**     | Count of oxygen or nitrogen atoms with lone pairs.   | Important for molecular recognition and binding.                                |
|      7 | **TPSA (Å²)**            | Topological Polar Surface Area.                      | Estimates a molecule’s ability to cross membranes (e.g., blood-brain barrier).  |
|      8 | **Rotatable Bonds**      | Number of flexible single bonds.                     | Higher values indicate more flexibility — affects oral bioavailability.         |
|      9 | **Formal Charge**        | Net ionic charge on the molecule.                    | Impacts solubility, reactivity, and membrane permeability.                      |
|     10 | **QED Score**            | Quantitative Estimate of Drug-likeness (0–1).        | A machine-learned score combining multiple rules — higher is better for drugs.  |
|     11 | **Solubility (LogS)**    | Aqueous solubility (in Log mol/L).                   | Critical ADME property — low solubility limits drug formulation.                |
|     12 | **Toxicity Probability** | Probability of being toxic to human cells (0–1).     | Prevents late-stage drug failure and ensures safety during drug discovery.      |

How are we predicitng?

| Property               | Type                   | Source   |
|------------------------|------------------------|----------|
| Canonical SMILES       | Rule-based             | RDKit    |
| Molecular Formula      | Rule-based             | RDKit    |
| Molecular Weight       | Rule-based             | RDKit    |
| LogP (Hydrophobicity)  | Rule-based             | RDKit    |
| H-Bond Donors          | Rule-based             | RDKit    |
| H-Bond Acceptors       | Rule-based             | RDKit    |
| TPSA (Polar Area)      | Rule-based             | RDKit    |
| Rotatable Bonds        | Rule-based             | RDKit    |
| Formal Charge          | Rule-based             | RDKit    |
| QED Score (Drug-likeness) | Heuristic ML-based  | RDKit    |
| Solubility (LogS)      | ML Regression Model    | ESOL + XGBoost |
| Toxicity Probability   | ML Classification Model| DeepChem + Tox21 |

---

## 🔧 Installation

1. Clone the repo:
   ```bash
   git clone https://github.com/therealsheero/SMILES-Property-Predictor.git
   cd smiles-property-predictor
````

2. Make sure your Python version is **3.10**. (Avoid 3.12+ for RDKit + DeepChem compatibility)

---

## 📂 Project Structure

```
.
├── Main_tool.ipynb           # Master script to run predictions
├── Solubility_Predictor.py     # ML model for solubility (LogS)
├── Toxicity_Predictor.py       # ML model for toxicity
├── models/
│   ├── logS_model_xgb_hybrid.pkl          # Pretrained solubility model
│   └── tox21_pretrained/       # DeepChem model for toxicity
└── Example_Sol/Tox.py          # Example script with sample SMILES
```

---

### ✅ Run the full property predictor:

```bash
python Main_tool.py
```

Enter any SMILES string when prompted (e.g. `CC(=O)Oc1ccccc1C(=O)O` for aspirin).

### 📄 Example output:

```
Canonical SMILES: CC(=O)Oc1ccccc1C(=O)O
Molecular Formula: C9H8O4
Molecular Weight: 180.159
LogP: 1.31
H-Bond Donors: 1
H-Bond Acceptors: 3
TPSA: 63.6 Å²
Rotatable Bonds: 2
Formal Charge: 0
QED Score: 0.55
Solubility (LogS): -2.29
Toxicity Probability: 0.399
```

---

## 📊 How the Models Work

* **Rule-based properties** are directly extracted using RDKit.
* **Solubility** is predicted using a **regression model** trained on the ESOL dataset with **XGBoost**.
* **Toxicity** is predicted using a **Graph Convolutional Network** from DeepChem trained on the **Tox21** dataset.

---

## 📚 References

* [RDKit Documentation](https://www.rdkit.org/)
* [ESOL Dataset](https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv)
* [Tox21 Toxicology Dataset](https://tripod.nih.gov/tox21/)
* [DeepChem](https://deepchem.io/)
* [QED Paper (Bickerton et al., 2012)](https://pubs.acs.org/doi/10.1021/jm300118s)

---


Let me know if you'd like me to add a `requirements.txt`, a Gradio or Streamlit UI version, or badge formatting at the top (`Built With RDKit · DeepChem · XGBoost`, etc).
```
