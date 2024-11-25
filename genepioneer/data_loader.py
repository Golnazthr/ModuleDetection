import os
import pandas as pd
from itertools import combinations
from collections import defaultdict

class DataLoader:
    def __init__(self, cancer_type):
        self.cancer_type = cancer_type
        self.TCGA_data_path = os.path.join("../GenesData/", self.cancer_type)
        self.IBM_data_path = os.path.join("../GenesData/IBP_GO_Terms.xlsx")


    def load_TCGA(self):
        genes_with_cases = defaultdict(set)

        genes_list_file_path = os.path.join(self.TCGA_data_path, f"{self.cancer_type}.tsv")
        genes_df = pd.read_csv(genes_list_file_path, sep='\t')
        genes_list = genes_df["Symbol"].tolist()
        # 👆 the list of 200 most significant genes for the chosen cancer type

        for gene in genes_list:
            gene_path = os.path.join(self.TCGA_data_path, gene)
            case_list_file_path = os.path.join(gene_path, f"{gene}.tsv")
            cases_df = pd.read_csv(case_list_file_path, sep='\t')
            cases_list = cases_df["Case ID"].tolist()

            # 👆 list of cases for each gene inside the list of genes for chosen cancer

            for case in cases_list:
                if "TCGA" in case: 
                    genes_with_cases[gene].add(case)
                    # 👆 list of cases that are related to TCGA will be appended to cases related to gene under process

        cases_with_genes = defaultdict(set)
        for gene, cases in genes_with_cases.items():
            for case in cases:
                cases_with_genes[case].add(gene)

        total_cases = len(cases_with_genes)
        
        return genes_with_cases, cases_with_genes, total_cases
    
    def load_IBM(self):
        processes_with_genes = defaultdict(set)

        processes_df = pd.read_excel(self.IBM_data_path)
        for index, row in processes_df.iterrows():
            process_name = row.iloc[0]
            genes_list = row.iloc[1:].dropna().tolist()
            processes_with_genes[process_name] = set(genes_list)

        genes_with_processes = defaultdict(set)

        for process, genes in processes_with_genes.items():
            for gene in genes:
                genes_with_processes[gene].add(process)               
  
        total_processes = len(processes_with_genes)
        return genes_with_processes, processes_with_genes, total_processes
        

