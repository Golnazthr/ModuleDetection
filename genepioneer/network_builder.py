import numpy as np
import networkx as nx
from collections import defaultdict
import csv

from .data_loader import DataLoader
from .network_analysis import NetworkAnalysis

class NetworkBuilder:
    def __init__(self, cancer_type):
        self.cancer_type = cancer_type
        self.graph = nx.Graph()
        data_loader = DataLoader(self.cancer_type)

        self.genes_with_cases ,self.cases_with_genes, self.total_cases = data_loader.load_TCGA()
        self.genes_with_processes, self.processes_with_genes, self.total_processes = data_loader.load_IBM()
                
        self.all_features = defaultdict(dict)

    def build_network(self):
        self.edge_adder()  
        return self.graph


    def edge_adder(self):
            frequent_genes = self.genes_with_cases.keys()
            for f_gene in frequent_genes:
                for process in self.genes_with_processes[f_gene]:
                    for gene in self.processes_with_genes[process]:
                        for gene_to_connect in self.processes_with_genes[process]:
                            if gene != gene_to_connect: 
                                if not self.graph.has_edge(gene, gene_to_connect):
                                    attributes = {
                                        'weight':  1
                                        }
                                    self.graph.add_edge(gene, gene_to_connect, **attributes)
                                else: 
                                    self.graph[gene][gene_to_connect]['weight'] += 1
    

    def weight_node(self, gi):
        total_weight = sum(data['weight'] for _, _, data in self.graph.edges(gi, data=True))
        return total_weight
    
    def weight_nodes(self, graph):
        return {node: self.weight_node(node) for node in graph.nodes()}

    def graph_entropy(self, weights):
        probabilities = weights / weights.sum()
        return -np.sum(probabilities * np.log(probabilities))

    def node_effect_on_entropy(self, node, entropy):
        copy_of_weights = self.node_weights_cache.copy()
        copy_of_weights.pop(node)
        
        for neighbor in self.graph.neighbors(node):
            if neighbor in copy_of_weights:
                edge_weight = self.graph[node][neighbor].get('weight', 0)
                copy_of_weights[neighbor] -= edge_weight
                copy_of_weights[neighbor] = max(copy_of_weights[neighbor], 1e-9)
        
        adjusted_weights = np.array(list(copy_of_weights.values()))
        new_entropy = self.graph_entropy(adjusted_weights)
        return abs(entropy - new_entropy)
    
    def calculate_all_features(self):
        closeness_centrality = nx.closeness_centrality(self.graph)
        print("closensess")
        
        betweenness_centrality = nx.betweenness_centrality(self.graph)
        print("betweenness")
        
        eigenvector_centrality = nx.eigenvector_centrality(self.graph)
        print("eigenvector")
        
        # Get weights for all nodes
        self.node_weights_cache = self.weight_nodes(self.graph)
        weights_array = np.array(list(self.node_weights_cache.values()))
        entropy = self.graph_entropy(weights_array)
        print("entropy")

        for node in self.graph.nodes():
            effect_on_entropy = self.node_effect_on_entropy(node, entropy)
            self.all_features[node] = {
                'weight': self.node_weights_cache[node],
                'closeness_centrality': closeness_centrality[node],
                'betweenness_centrality': betweenness_centrality[node],
                'eigenvector_centrality': eigenvector_centrality[node],
                'effect_on_entropy': effect_on_entropy,
            }
            self.graph.nodes[node].update({
                'weight': self.node_weights_cache[node],
                'closeness_centrality': closeness_centrality[node],
                'betweenness_centrality': betweenness_centrality[node],
                'eigenvector_centrality': eigenvector_centrality[node],
                'effect_on_entropy': effect_on_entropy,
                'graph_entropy': entropy
            })
        
        print("goes to ls")
        self.add_LS_to_network()
        
        self.all_features['graph_entropy'] = entropy
        
        return self.all_features
    
    def save_features_to_csv(self, all_features, filename):
        nx.write_gml(self.graph, f"{filename}.gml")
        with open(f"{filename}.csv", 'w', newline='') as csvfile:
            fieldnames = ['node', 'weight', 'closeness_centrality', 'betweenness_centrality',
                        'eigenvector_centrality', 'effect_on_entropy', "ls_score"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for node, features in all_features.items():
                if node != 'graph_entropy' and node != "ls_for_features":
                    features['node'] = node  # Add node ID to the row
                    writer.writerow(features)
                    
            writer.writerow({'node': 'graph_entropy', 'weight': all_features['graph_entropy']})
            writer.writerow({'node': 'ls_for_features', 'weight': all_features['ls_for_features']})
                
            
    def add_LS_to_network(self):
        network_analysis = NetworkAnalysis(self.cancer_type, self.all_features);

        nodes, features = network_analysis.compute_laplacian_scores()
        for node in self.graph.nodes():
                        
            self.all_features[node]["ls_score"] = nodes.loc[node].LaplacianScore
            
            self.graph.nodes[node].update({
                'ls_score': nodes.loc[node].LaplacianScore,
                'ls_for_features': np.array_str(features),
            })
            
        self.all_features["ls_for_features"] = np.array_str(features)
        return self.all_features
    