
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
   git clone https://github.com/yourusername/smiles-property-predictor.git
   cd smiles-property-predictor
````

2. Install required dependencies:

   ```bash
   pip install -r requirements.txt
   ```

3. Make sure your Python version is **3.10**. (Avoid 3.12+ for RDKit + DeepChem compatibility)

---

## 📂 Project Structure

```
.
├── main_predictor.py           # Master script to run predictions
├── Solubility_Predictor.py     # ML model for solubility (LogS)
├── Toxicity_Predictor.py       # ML model for toxicity
├── models/
│   ├── logS_model.pkl          # Pretrained solubility model
│   └── tox21_pretrained/       # DeepChem model for toxicity
└── test_run.py                 # Example script with sample SMILES
```

---

## ▶️ How to Use

### ✅ Run the full property predictor:

```bash
python main_predictor.py
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

## 🙌 Future Add-ons

* Web dashboard (Streamlit/Gradio)
* Batch SMILES upload and export
* More ML-based properties (e.g., bioactivity, synthetic accessibility)

---

## 🧠 Author

**Sheero** — Exploring AI in Drug Discovery and Bioinformatics 🚀
Let’s build BioMedGPT together!

---

## 📜 License

MIT License. Free to use and contribute!

```

---

Let me know if you'd like me to add a `requirements.txt`, a Gradio or Streamlit UI version, or badge formatting at the top (`Built With RDKit · DeepChem · XGBoost`, etc).
```
