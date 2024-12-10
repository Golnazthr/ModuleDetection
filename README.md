# GenePioneer: A Comprehensive Python Package for Identification of Essential Genes and Modules in Cancer


This repository contains two folders: **GenePioneer**, a Python package for gene analysis, and **Genes-Data**, a collection of data gathered using a web scraper.


## Description

The GenePioneer was developed as a fast and straightforward way to integrate gene ranking and module detection into a practical, Python-based tool for cancer researchers. It requires minimal input, delivers clear output,
and can be run within a Python environment, making it highly user-friendly and accessible to non expert programmers while supporting large-scale dataset analysis. By evaluating gene importance and identifying gene interactions within cancer
networks, GenePioneer provides critical insights into the genetic drivers of cancer. Key features include ranking genes by their network significance and identifying the modules they belong to, which helps explore cancer-related pathways and aids in developing precise therapies. GenePioneer’s user-centric design ensures that researchers of all skill levels can make use of its capabilities. By combining comprehensive data integration, advanced networkbased analysis, and statistical rigor, GenePioneer stands as a versatile and impactful resource for cancer research across multiple cancer types.

![Image description](Workflow.pdf)

### Questions about the implementation:

Amirhossein Haerianardakani, haerian.amirhossein[at]gmail.com


### Installation

Install the package from PyPI:

```bash
pip install genepioneer
```

### Usage

Follow these steps to use the Gene Pioneer package for gene analysis:

1. Create a Python File: Create a Python file in your local environment, for example gene_analyze.py.

2. Import the Gene Pioneer Class: Import the main class designed for gene analysis:

```bash
from genepioneer import GeneAnalysis`
```

3. Prepare a Gene List: Create a `.txt` file containing a list of genes in the `OFFICIAL_GENE_SYMBOL` format. Each line should contain one gene name. For example:

```
BRCA1
TP53
PTEN
```

4. Initialize the GeneAnalysis Class: Pass the cancer type and the path to the gene list file to the constructor.

```bash
gene_analysis = GeneAnalysis("Ovary", "./benchmark-data/gene_list.txt")
```

- Supported Cancer Types:

`"Adrenal", "Bladder", "Brain", "Cervix", "Colon", "Corpus uteri", "Kidney", "Liver", "Ovary", "Prostate", "Skin", "Thyroid"`

- Example Input:
  `Ovary` is the cancer type.
  `gene_list.txt` is the gene list file.

5. Run Gene Analysis:

```bash
gene_analysis.analyze_genes()
```

After running this, an `output.json` file will be generated in your current working directory.

### Code Overview

Here is a summary of the main files in the genepioneer package:

`__init__.py`: Initializes the Gene Pioneer package.
`data_loader.py`: Handles data loading operations.
`evaluation.py`: Contains evaluation metrics and methods.
`gene_analysis.py`: Implements the main GeneAnalysis class for gene processing.
`network_analysis.py`: Provides network-based analysis for gene interactions.
`network_builder.py`: Builds networks for gene interaction analysis.
