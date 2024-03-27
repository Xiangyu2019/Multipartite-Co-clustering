# Network Analysis and Drug Repurposing Framework

This repository contains a collection of Python scripts designed for the analysis of biological networks, drug repurposing, and exploration of drug-target interactions within the context of protein-protein interactions (PPI), gene interactions (GI), and molecular interactions (MI). It leverages data parsing, network smoothing, topological analysis, and matrix computation techniques to uncover insights into drug similarities and potential repurposing opportunities.

## Description

The framework is structured around several key scripts, each serving a specific purpose in the process of discovering protein interactions and repurposing drugs:

## Data Processing and Matrix Creation

- `scripts/Data/computeDrugSimilarities.py`: Computes and filters Tanimoto similarities between approved and experimental drugs based on their SMILES representation, optionally saving the results.

- `scripts/Data/creatingMatrices.py`: Handles the creation, augmentation, and initialization of matrices derived from biological data, focusing on interactions between SARS-CoV-2 proteins, human proteins, and drugs, utilizing BioGRID and DrugBank datasets.

- `scripts/Data/creatingNetworksUtils.py`: Transforms and manages network data, including converting edge lists to adjacency and matrices, filtering repeated interactions, and processing raw BioGRID data into network form, with options for undirected graphs and species-specific filtering.

- `scripts/Data/parsingDrugBank.py`: Converts the data into structured CSV formats for different aspects of drug information, including groups, SMILES representations, targets, and categories, with options to save the extracted data.

- `scripts/Data/smoothingBipartiteNetwork.py`: Spreads the influence of nodes in a bipartite matrix across a network to achieve a smooth distribution of values. It provides functionality to normalize a network, apply the iterative propagation process with a stop criterion based on the Frobenius norm, and save the results for different values of the alpha parameter, which controls the balance between propagation and original values.

- `scripts/Data/constructing the Metabolic Interaction Network/create_KEGG_metabolic_modified.sh`: Automates the process of fetching metabolic pathway data from the KEGG database for Homo sapiens (human), extracting specific pathways related to metabolism, and then downloading detailed information for each of those pathways.

## Data Fusion Framework

- `scripts/DataFusionFramework/choosingParameters.py`: Analyzes network structures across various parameter configurations, allowing the evaluation of clustering consistency across multiple runs.

- `scripts/DataFusionFramework/framework.py`: Implements the data fusion framework aimed at integrating different types of biological data represented as matrices. The function takes as input the relational matrices between different entities (e.g., genes, diseases, drugs), the Laplacian matrices capturing the relationships within entities, and initial guesses for the factor matrices.

## Results Analysis and Visualization

- `scripts/ResultsAnalysis/PRandROCcurves.py`: Contains functions for the evaluation and application of drug-target interaction (DTI) predictions using matrix factorization techniques.

- `scripts/ResultsAnalysis/comparisonConstituentMatrices.py`: Takes a pandas DataFrame representing an edge list, where each row contains two nodes forming an edge. It generates a new list of edges by concatenating node pairs into strings, ensuring a standardized format regardless of the node order.

- `scripts/ResultsAnalysis/enrichmentAnalysis.py`: A collection of functions designed for cluster analysis, annotation loading, and enrichment analysis within biological datasets, particularly focusing on genes and drugs.

- `scripts/ResultsAnalysis/topologicalAnalysisNetwork.py`: Processes biological networks to identify gene sets based on connectivity and expression data, calculates various centrality measures to analyze network structure, and saves results for further analysis.

## Plotting and Visualization

- `scripts/plots.py`: Visualizes data analysis results, including plotting enrichment percentages for genes and drugs, association scores, precision-recall, and ROC curves for DTIs.

## Jupyter Notebooks

- `Comparison constituent networks (PPI, GI MI) with MIN.ipynb`: Processes and analyzes biological networks represented as edge lists, focusing on Protein-Protein Interactions (PPI), Genetic Interactions (GI), Metabolic Interactions (MI), and a combined Multilayer Interaction Network (MIN).

- `Enrichment Analysis of the Clusters of Genes and Drugs - MIN and PPI.ipynb`: Performs enrichment analyses on drug and gene clusters derived from a data fusion framework, calculating the percentage of enriched categories, clusters, and entities.

- `Executing data fusion framework - Parametrization Tunning - Dispersion Coefficient.ipynb`: Computes the coefficient for various parameter configurations in a data fusion framework, assessing clustering quality from viral proteins, genes, and drugs data to identify optimal clustering parameters.

- `Executing data fusion framework - Parametrization Tunning.ipynb`: Implements the data fusion framework to integrate matrices representing interactions between different biological entities (viral proteins, genes, and drugs). It explores various parameter settings (like the number of clusters for each entity type) across multiple runs, initializing matrices randomly based on a set of predefined seeds to ensure reproducibility.

- `Executing data fusion framework.ipynb`: Executes the data fusion process, focusing on integrating biological interaction data into a unified framework. It includes conditional paths for creating or loading matrices based on whether the PPI network is treated as a comprehensive human interactome. After determining or loading the necessary matrices (including relational matrices and Laplacians for interaction networks), the script proceeds with a specified parametrization to initialize matrices for the fusion process. it integrates these matrices, aiming to reveal complex interrelations between viral proteins, human genes, and drugs.

- `Matrix Completion Property of the DTI matrix - Choosing the threshold - PR and ROC - MIN and PPI.ipynb`: A comprehensive analysis of reconstructed drug-target interaction (DTI) matrices for both a Metabolic Interaction Network (MIN) and a Protein-Protein Interaction (PPI) network treated as the human interactome. It begins by loading the reconstructed DTI matrices and their original counterparts, then calculates precision-recall (PR) and receiver operating characteristic (ROC) curves to evaluate the models' performance.

- `Topological Analysis Network.ipynb`: Analyzes viral infection and differentially expressed genes within protein-protein and metabolic interaction networks. It creates gene sets, assesses overlap, tests for statistical significance, examines degree distributions, and compares centrality measures.

- `Executing data fusion framework - Parametrization Tunning.ipynb`: The process of creating, tuning, and saving matrices for a data fusion framework.


## Getting Started

### Dependencies

- Python 3.8
- NumPy, Pandas, SciPy, NetworkX
- Requirements.txt
  
### Installation

To set up the project, start by cloning the repository to your local machine. After cloning, you will need to install the required Python packages to ensure the project runs smoothly. Additionally, a shell script that assists in constructing the network is available and can be found at `scripts/Data/constructing the Metabolic Interaction Network/create_KEGG_metabolic_modified.sh`.


## Data

The necessary data for executing the framework's scripts is readily accessible in the Hugging Face repository. You can find all required datasets at (https://huggingface.co/datasets/Xiangyuden/Network-Fusion). This ensures integration and usage within the framework.

## Contributing

We welcome contributions to improve the framework and extend its functionalities. Please feel free to fork the repository, make your changes, and submit a pull request.



