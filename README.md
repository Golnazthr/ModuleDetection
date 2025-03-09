# GenePioneer: A Comprehensive Python Package for Identification of Essential Genes and Modules in Cancer


## Description

<img src="Workflow.png" width="40%" align="right" />

The GenePioneer was developed as a fast and straightforward way to integrate gene ranking and module detection into a practical, Python-based tool for cancer researchers. It requires minimal input, delivers clear output, and can be run within a Python environment, making it highly user-friendly and accessible to non expert programmers while supporting large-scale dataset analysis. By evaluating gene importance and identifying gene interactions within cancer networks, GenePioneer provides critical insights into the genetic drivers of cancer. Key features include ranking genes by their network significance and identifying the modules they belong to, which helps explore cancer-related pathways and aids in developing precise therapies.  GenePioneer’s user-centric design ensures that researchers of all skill levels can make use of its capabilities. By combining comprehensive data integration, advanced networkbased analysis, and statistical rigor, GenePioneer stands as a versatile and impactful resource for cancer research across multiple cancer types.

### Features

- **Gene Ranking**: Determines gene importance based on network significance.
- **Module Detection**: Identifies gene clusters within cancer pathways.
- **Statistical Analysis**: Evaluates detected gene modules and their association with known pathways.
- **User-Friendly API**: Allows researchers of all skill levels to analyze cancer-related genetic data efficiently.
- **Full Reproducibility**: Option to either use precomputed data or regenerate all components from raw datasets.

## Installation

GenePioneer is available via PyPI. Install it using:

```bash
pip install genepioneer
```

## Two Usage Modes

You can either:

1. **Use Preprocessed Data**: Run analysis using prebuilt datasets in `Data/cancer-gene-data/` and `Data/module-data/`.
2. **Reproduce Everything**: Build networks, generate rankings, and detect modules from raw cancer data stored in `GenesData/`.

---

## **Option 1: Using Preprocessed Data**

This is the simplest approach. You only need to provide a list of genes and specify the cancer type.

### **Step 1: Prepare a Gene List**

Create a `.txt` file containing one gene name per line in the **OFFICIAL_GENE_SYMBOL** format.

#### Example (`gene_list.txt`):

```
BRCA1
TP53
PTEN
```

### **Step 2: Run Gene Analysis**

```python
from genepioneer import GeneAnalysis

gene_analysis = GeneAnalysis("Ovary", "./Data/benchmark-data/gene_list.txt")
gene_analysis.analyze_genes()
```

### **Step 3: Output**

This will generate an `output.json` file with:

- **Gene Rankings**: Sorted based on importance in the network.
- **Modules**: Groups of genes functionally related in cancer.
- **Statistical Significance**: Evaluation of identified modules.

#### **Supported Cancer Types**

```
"Adrenal", "Bladder", "Brain", "Cervix", "Colon", "Corpus uteri", "Kidney", "Liver", "Ovary", "Prostate", "Skin", "Thyroid"
```

---

## **Option 2: Reproducing Everything (Building Data from Scratch)**

If you want full control over data generation, follow these steps to build your own cancer-specific datasets.

### **Step 1: Add Required Data**

You need:

- **Raw TCGA Cancer Data (`GenesData/`)**: Cancer-specific gene expression data.
- **IBM Gene Ontology (`GenesData/IBP_GO_Terms.xlsx`)**: Gene-to-biological process mappings.

#### **Example Directory Structure**

```
GenesData/
│-- IBP_GO_Terms.xlsx
│-- Adrenal/
│   │-- ABL1/
│   │   │-- ABL1.tsv
```

### **Step 2: Build Network and Compute Features**

```python
from genepioneer import NetworkBuilder

network_builder = NetworkBuilder("Adrenal", "./GenesData")
graph = network_builder.build_network()
features = network_builder.calculate_all_features()
network_builder.save_features_to_csv(features, "./Data/cancer-gene-data/Adrenal")
```

This step:

- Builds a **gene interaction network**.
- Computes **network-based features** (e.g., centrality, entropy, Laplacian scores).
- Saves features to a CSV file.

### **Step 3: Detect Modules**

```python
from genepioneer import NetworkAnalysis

network_analysis = NetworkAnalysis("Adrenal", features)
modules = network_analysis.module_detection()
```

- Identifies **gene modules** based on connectivity and functional relevance.
- Saves results as `Data/module-data/Adrenal.json`.

### **Step 4: Run Full Gene Analysis**

Once networks and modules are generated, you can proceed with standard gene analysis:

```python
from genepioneer import GeneAnalysis

gene_analysis = GeneAnalysis("Adrenal", "./Data/benchmark-data/gene_list.txt", 
                             cancer_gene_path="./Data/cancer-gene-data", 
                             module_data_path="./Data/module-data")
gene_analysis.analyze_genes()
```

---

## **Dataset Structure and Format**

### **1. TCGA Data (`GenesData/*/*.tsv`)**

Contains gene expression and associated cases.

#### Example (`GenesData/Adrenal/ABL1/ABL1.tsv`):

```
Case ID	Expression
TCGA-01	2.5
TCGA-02	1.8
```

### **2. IBM Gene Ontology (`GenesData/IBP_GO_Terms.xlsx`)**

Links genes to biological processes.

| Process | Gene1 | Gene2 |
|---------|------|------|
| Cell Cycle | BRCA1 | TP53 |

### **3. Network Features (`Data/cancer-gene-data/*.csv`)**

Stores computed network importance scores.

#### Example (`Data/cancer-gene-data/Adrenal_network_features.csv`):

```
node,ls_score
ABL1,0.85
TP53,0.92
```

### **4. Module Data (`Data/module-data/*.json`)**

Contains detected gene modules.

#### Example (`Data/module-data/Adrenal.json`):

```json
{
  "module_1": [
    ["ABL1", "TP53"],
    3.5,
    1.2
  ]
}
```

## **Reproducibility Steps**

1. **Clone the Repository**

```bash
git clone https://github.com/yourusername/GenePioneer.git
cd GenePioneer
```

2. **Install Dependencies**

```bash
pip install -r requirements.txt
```

3. **Add or Generate Data**

- **Use prebuilt data** (Option 1), or
- **Generate data from raw sources** (Option 2).

4. **Run Gene Analysis**

```bash
python -m genepioneer.gene_analysis "Adrenal" "./Data/benchmark-data/gene_list.txt"
```

5. **Verify Output**

- `output.json` contains ranked genes and detected modules.


## Questions about the implementation:

Amirhossein Haerianardakani, haerian.amirhossein[at]gmail.com


If you encounter a bug, experience a failed function, or have a feature request, please open an issue in the GitHub or contact Amirhossein.

## License

This project is licensed under the MIT License - [MIT License](https://opensource.org/licenses/MIT)
