#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd


# In[ ]:


get_ipython().run_line_magic('run', '../scripts/ResultsAnalysis/topologicalAnalysisNetwork.py')
get_ipython().run_line_magic('run', '../scripts/Data/creatingNetworksUtils.py')
get_ipython().run_line_magic('run', '../scripts/ResultsAnalysis/comparisonConstituentMatrices.py')
get_ipython().run_line_magic('run', '../scripts/plots.py')


# ### Load datasets

# In[ ]:


edgeListPPI = pd.read_csv('../Data/Networks/BioGRID_PPI-HS_edgeList.csv', index_col=0)
dfPPI = createAdjacencyMatrixFromEdgeList(edgeListPPI, undirected=True)


# In[ ]:


edgeListGI = pd.read_csv('../Data/Networks/BioGRID_GI-HS_edgeList.csv', index_col=0)
dfGI = createAdjacencyMatrixFromEdgeList(edgeListGI, undirected=True)


# In[ ]:


edgeListMI = pd.read_csv('../Data/Networks/KEGG_metabolic_network_hgnc.csv', index_col=0)
dfMI = createAdjacencyMatrixFromEdgeList(edgeListMI, undirected=True)


# In[ ]:


edgeListMIN = pd.read_csv('../Data/Networks/MIN_edgeList.csv')
dfMIN = createAdjacencyMatrixFromEdgeList(edgeListMIN, undirected=True)


# ### Gene ovelap between the consituent networks

# In[ ]:


genesPPI = dfPPI.columns
genesGI = dfGI.columns
genesMI = dfMI.columns


# In[ ]:


fileName = '../Plots/overlapGenes_consituentMatrices.png'
plotThreeOverlapSets(genesPPI, genesGI, genesMI, fileName)


# ### Edge Overlap

# In[ ]:


edgesPPI = convertEdgeList(edgeListPPI)
edgesGI = convertEdgeList(edgeListGI)
edgesMI = convertEdgeList(edgeListMI)


# In[ ]:


fileName = '../Plots/overlapEdges_consituentMatrices.png'
plotThreeOverlapSets(edgesPPI, edgesGI, edgesMI, fileName)


# ### Basic network centrality

# In[ ]:


centralitiesPPI = centralityMeasures(dfPPI, verbose=True)
centralitiesGI = centralityMeasures(dfGI, verbose=True)
centralitiesMI = centralityMeasures(dfMI, verbose=True)
centralitiesMIN = centralityMeasures(dfMIN, verbose=True)


# In[ ]:


dfCentralities = pd.DataFrame([centralitiesPPI, centralitiesGI, centralitiesMI, centralitiesMIN], index=['PPI', 'GI', 'MI', 'MIN'], columns=['Average Degree', 'Eigenvector Centrality', 'Clustering Coefficient', 'Betweeness Centrality', 'Closeness Centrality'])
dfCentralities.to_csv('../Data/Results/ConstituentMatrices_centarlitiesComparison.csv')


# ## Compare GDV of the three constituent networks

# #### Load networks

# In[ ]:


edgesListMIN = pd.read_csv('../Data/Networks/MIN_edgeList.csv', header=None)
MIN = createAdjacencyMatrixFromEdgeList(edgesListMIN)


# In[ ]:


PPI = loadData('../Data/Networks/BioGRID_PPI-HS_edgeList.csv')
GI = loadData('../Data/Networks/BioGRID_GI-HS_edgeList.csv')
MI = loadData('../Data/Networks/KEGG_metabolic_network_hgnc.csv')


# #### Create orca files

# In[ ]:


label_mappingPPI, reverse_mappingPPI = createOrcaFile(PPI, "PPI_names_edgelist")
label_mappingGI, reverse_mappingGI = createOrcaFile(GI, "GI_names_edgelist")
label_mappingMI, reverse_mappingMI = createOrcaFile(MI, "MI_names_edgelist")
label_mappingMIN, reverse_mappingMIN = createOrcaFile(MIN, "MIN_names_edgelist")


# #### Plot

# In[ ]:


saveFile = '../Plots/GDV_constituentMatrices.png'
outFilePath = '../Data/OrcaFormatFiles/'
plotGDVDifferentNetworks(["PPI_names_edgelist", "GI_names_edgelist", "MI_names_edgelist", "MIN_names_edgelist"], outFilePath, [reverse_mappingPPI, reverse_mappingGI, reverse_mappingMI, reverse_mappingMIN], ['PPI', 'GI', 'MI', 'MIN'], saveFile)

