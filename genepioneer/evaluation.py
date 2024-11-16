import pandas as pd
import glob
import os

import concurrent.futures


import numpy as np

import json
from gprofiler import GProfiler

from sklearn.metrics import precision_score, recall_score, roc_curve, auc
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


from genepioneer import DataLoader

class Evaluation:
    def __init__(self):
        self.gp = GProfiler(return_dataframe=True)
        
        self.benchmark_genes = self.read_benchmark_genes("../Data/benchmark-data")
        self.network_genes = self.read_network_genes("../Data/cancer-gene-data")
        self.modules = self.read_modules("../Data/module-data")
        self.result = self.eval(self.network_genes, self.benchmark_genes)
        self.module_results = self.evaluate_modules(self.modules)
    
    def read_modules(self, benchmark_folder):
        benchmark_genes = {}
        for filepath in glob.glob(os.path.join(benchmark_folder, '*.json')):
            module_name = os.path.basename(filepath).split('.')[0]
            with open(filepath, 'r') as file:
                modules = json.load(file)
            benchmark_genes[module_name] = modules
        return benchmark_genes
                
    def read_benchmark_genes(self, benchmark_folder):
        benchmark_genes = {}
        for filepath in glob.glob(os.path.join(benchmark_folder, '*.txt')):
            benchmark_name = os.path.basename(filepath).split('.')[0]
            with open(filepath, 'r') as file:
                genes = set(file.read().strip().split('\n'))
            benchmark_genes[benchmark_name] = genes
        return benchmark_genes

    def read_network_genes(self, network_folder):
        network_genes = {}
        for filepath in glob.glob(os.path.join(network_folder, '*.csv')):
            cancer_type = os.path.basename(filepath).split('.')[0]
            df = pd.read_csv(filepath)
            df = df.sort_values(by='ls_score', ascending=False)
            genes = df['node'].tolist()
            network_genes[cancer_type] = genes
        return network_genes
    def read_mutated_genes(self, network_folder):
        network_genes = {}
        for filepath in glob.glob(os.path.join(network_folder, '*.csv')):
            cancer_type = os.path.basename(filepath).split('.')[0]
            cancer_type = cancer_type.replace('_network_features', '')
            data_loader = DataLoader(cancer_type)
            genes_with_cases, cases_with_genes, total_cases = data_loader.load_TCGA()
            genes = list(genes_with_cases.keys())
            network_genes[cancer_type] = genes
        return network_genes

    def get_top_n_genes(self, genes, n):
        return set(genes[:n])

    def calculate_metrics(self, predicted_genes, benchmark_genes):
        TP = len(predicted_genes.intersection(benchmark_genes))
        B_size = len(benchmark_genes)
                
        return TP, B_size

    def eval(self, network_genes, benchmark_genes):
        results = {}
        for cancer_type, genes in network_genes.items():
            results[cancer_type] = {}
            total_genes = set(genes)
            for benchmark_name, benchmark in benchmark_genes.items():
                top_n_genes = self.get_top_n_genes(genes, len(benchmark))
                metrics = self.calculate_metrics(top_n_genes, benchmark)
                auc_roc, precision, recall = self.calculate_auc_roc(top_n_genes, benchmark, total_genes)
                driver_ranks, other_ranks = self.evaluate_driver_genes(genes, benchmark)
                results[cancer_type][benchmark_name] = {
                    'metrics': metrics,
                    'auc_roc': auc_roc,
                    'precision': precision,
                    'recall': recall,
                    'driver_ranks': driver_ranks,
                    'other_ranks': other_ranks
                }
        return results
    
    def evaluate_driver_genes(self, genes, benchmark):
        driver_ranks = [genes.index(gene) for gene in benchmark if gene in genes]
        other_ranks = [rank for rank, gene in enumerate(genes) if gene not in benchmark]
        return driver_ranks, other_ranks

    def calculate_auc_roc(self, predicted_genes, benchmark_genes, total_genes):
        y_true = [1 if gene in benchmark_genes else 0 for gene in total_genes]
        y_scores = [1 if gene in predicted_genes else 0 for gene in total_genes]
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        precision = precision_score(y_true, y_scores)
        recall = recall_score(y_true, y_scores)
        auc_roc = auc(fpr, tpr)
        return auc_roc, precision, recall
    
    def print_driver_ranking_stats(self, cancer_type, benchmark_name, driver_ranks, other_ranks):
        if driver_ranks:
            mean_rank = np.mean(driver_ranks)
            median_rank = np.median(driver_ranks)
            rank_percentiles = np.percentile(driver_ranks, [25, 50, 75])
            mannwhitney_p = mannwhitneyu(driver_ranks, other_ranks, alternative='less').pvalue
            print(f'Cancer Type: {cancer_type} - Benchmark: {benchmark_name}')
            print(f'  Median Rank: {median_rank}')
            print(f'  Rank Percentiles (25th, 50th, 75th): {rank_percentiles}')
            print(f'  Mann-Whitney U Test p-value: {mannwhitney_p}')
        else:
            print(f'Cancer Type: {cancer_type} - Benchmark: {benchmark_name}')
            print('  No driver genes found in the ranked list.')
    
    def print_result(self):
        for cancer_type, benchmarks in self.result.items():
            print(f'Cancer Type: {cancer_type}')
            for benchmark_name, result in benchmarks.items():
                metrics = result['metrics']
                TP, B_size = metrics
                auc_roc = result['auc_roc']
                precision = result['precision']
                recall = result['recall']
                driver_ranks = result['driver_ranks']
                other_ranks = result['other_ranks']
                print(f'Benchmark: {benchmark_name}')
                print(f'TP: {TP}, size of benchmark: {B_size}')
                print(f'AUC-ROC: {auc_roc:.3f}')
                print(f'precision: {precision:.3f}')
                print(f'recall: {recall:.3f}')
                self.print_driver_ranking_stats(cancer_type, benchmark_name, driver_ranks, other_ranks)
                
    def evaluate_modules(self, modules):
        
        pathways_of_interest = {
            "KEGG:04068",
            "KEGG:04310",
            "KEGG:04010",
            "KEGG:04115",
            "KEGG:04915",
            "KEGG:04014",
            "KEGG:04012",
            "KEGG:04150",
            "KEGG:05200",
            "KEGG:04151",
            "KEGG:04370"
        }
        results = {}        
        def evaluate_single_module(module):
            genes, score1, score2 = module
            enrichment_results = self.gp.profile(organism='hsapiens', query=genes)
            enrichment_results.to_csv('out.csv')
            significant_pathways = enrichment_results[
                (enrichment_results['p_value'] <= 0.05) & 
                (enrichment_results['native'].isin(pathways_of_interest))
            ]
            if len(significant_pathways) >= 2:
                return {
                    'module_genes': genes,
                    'score1': score1,
                    'score2': score2,
                    'significant_pathways': significant_pathways[['name', 'p_value']].to_dict(orient='records')
                }
            
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_to_cancer_type = {
                executor.submit(evaluate_single_module, module): (cancer_type, module)
                for cancer_type, module_list in modules.items() for module in module_list
            }
            
            count = 0
            for future in concurrent.futures.as_completed(future_to_cancer_type):
                cancer_type, module = future_to_cancer_type[future]
                if cancer_type not in results:
                    results[cancer_type] = []
                try:
                    evaluation = future.result()
                    if evaluation: 
                        results[cancer_type].append(evaluation)
                    count += 1
                    print(count)
                except Exception as e:
                    print(f"Module evaluation failed for {cancer_type}: {e}")
    
        return results
    
    def print_module_evaluation(self): 
        with open('evaluated_modules_result.json', 'w') as f:
            json.dump(self.module_results, f, indent=2)
