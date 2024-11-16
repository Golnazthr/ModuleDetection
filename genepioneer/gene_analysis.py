import json
import pandas as pd

import networkx as nx
import pkg_resources

class GeneAnalysis:

    def __init__(self, cancer_type, input_path):
        self.cancer_type = cancer_type
        
        self.genes = self.load_genes()
        self.modules = self.load_modules()
        self.rankings = self.rank_genes_by_ls_score()
        self.test_genes = self.load_test_genes(input_path)
                

    def load_genes(self):
        data_path = pkg_resources.resource_filename(__name__, f"Data/cancer-gene-data/{self.cancer_type}_network_features.csv")
        df = pd.read_csv(data_path)

        gene_data = df.to_dict(orient='records')
    
        return gene_data
    
    def load_modules(self):
        data_path = pkg_resources.resource_filename(__name__, f"Data/module-data/{self.cancer_type}.json")
        with open(data_path, 'r') as file:
            data = json.load(file)
        
        modules = {}
        for module_index, (gene_list, score, additional_value) in enumerate(data):
            for gene in gene_list:
                if gene not in modules:
                    modules[gene] = []
                modules[gene].append({'genes': gene_list, 'score': score, 'additional_value': additional_value})
        
        return modules
    
    def load_test_genes(self, input_path):
        # Read genes from the specified input file
        with open(input_path, 'r') as file:
            genes = [line.strip() for line in file if line.strip()]
        return genes
    
    def rank_genes_by_ls_score(self):
    # Sort the gene data by 'ls_score' in descending order and assign ranks
        ranked_genes = {gene['node']: rank for rank, gene in enumerate(sorted(self.genes, key=lambda x: x.get('ls_score', 0), reverse=True))}

        return ranked_genes
    
    def find_gene_modules(self, gene):
        return self.modules.get(gene, [])
    
    def analyze_genes(self):
        results = {}
        for gene in self.test_genes:
            gene_info = {}
            gene_info['rank'] = self.rankings.get(gene, 'Not in network')
            gene_info['modules'] = self.find_gene_modules(gene)
            results[gene] = gene_info
        
        with open("./output.json", 'w') as f:
            json.dump(results, f, indent=4)