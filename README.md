# Network Analysis and Drug Repurposing Framework

This repository contains a collection of Python scripts designed for the analysis of biological networks, drug repurposing, and exploration of drug-target interactions within the context of protein-protein interactions (PPI), gene interactions (GI), and molecular interactions (MI). It leverages data parsing, network smoothing, topological analysis, and matrix computation techniques to uncover insights into drug similarities and potential repurposing opportunities.

## Description

The framework is structured around several key scripts, each serving a specific purpose in the process of discovering protein interactions and repurposing drugs:

- **Data Parsing and Preparation**: Utilize `parsingDrugBank.py`, `creatingMatrices.py`, and `creatingNetworksUtils.py` for initial data handling and matrix generation.
- **Network Smoothing and Analysis**: Apply `smoothingBipartiteNetwork.py` for network smoothing, and `Topological Analysis Network.py` for understanding network structures. Additionally, `Comparison constituent networks (PPI, GI MI) with MIN.py` aids in comparing network types.
- **Matrix and Computational Analysis**: Leverage `computeDrugSimilarities.py` for drug similarity computations, essential for drug repurposing insights.
- **Enrichment Analysis**: Use `Enrichment Analysis of the Clusters of Genes and Drugs - MIN and PPI.py` along with `Enrichment Analysis Gene Sets in MIN.py` for deep dives into gene and drug clusters.
- **Data Fusion and Parameter Optimization**: Engage with `Executing data fusion framework.py` and its related parameter tuning scripts to refine and optimize the predictive models.

## Getting Started

### Dependencies

- Python 3.8
- NumPy, Pandas, SciPy, NetworkX
  
### Installation

Clone the repository and install the required Python packages.


## Data

The data essential for executing the framework's scripts can be found in the repo

## Contributing

We welcome contributions to improve the framework and extend its functionalities. Please feel free to fork the repository, make your changes, and submit a pull request.



