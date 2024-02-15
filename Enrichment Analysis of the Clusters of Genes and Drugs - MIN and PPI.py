#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('run', '../scripts/ResultsAnalysis/enrichementAnalysis.py')
get_ipython().run_line_magic('run', '../scripts/plots.py')


# In[ ]:


decompositionPPI = False
if decompositionPPI:
    pathFile = '../Data/Results/PPIasHumanInteractome/'
    genes = pd.read_csv('../Data/Names/Genes_PPI.csv')['Genes'].values
    saveFig = '../Plots/ClusteringEnrichmentsPPI.png'
else:
    pathFile = '../Data/Results/'
    genes = pd.read_csv('../Data/Names/Genes.csv')['Genes'].values
    saveFig = '../Plots/ClusteringEnrichments.png'


# ### Drugs

# In[ ]:


#Load drugs
drugs = pd.read_csv('../Data/Names/Drugs.csv')['Drugs'].values

#Extract clusters
clusters = extractClusters(drugs, pathFile+'G3_SVD.npy')

#Load annotations data
categoriesToDrugs = loadDrugAnnotations(drugs, verbose=True)

#Compute enrichment Analysis
categoriesEnriched, clustersEnriched, drugsEnriched = enrichmentAnalysis(clusters, categoriesToDrugs)

dfDrugsEnrichments = pd.DataFrame([categoriesEnriched, clustersEnriched, drugsEnriched], index=['Clusters of Drugs', 'Drug Categories', 'Drugs'], columns=['DC']).transpose()
dfDrugsEnrichments


# ### Genes

# In[ ]:


#Extract clusters
clusters = extractClusters(genes, pathFile+'G2_SVD.npy')
    
GOEnrichments = []
aspectsLabel = []
for aspect in ['BiologicalProcess', 'MolecularFunction', 'CellularComponent']:
    #Load annotations data
    annotationsToGenes = loadGeneAnnotations(aspect, genes, verbose=False)

    #Compute enrichment Analysis
    annotationsEnriched, clustersEnriched, genesEnriched = enrichmentAnalysis(clusters, annotationsToGenes)

    GOEnrichments.append([annotationsEnriched, clustersEnriched, genesEnriched])
    
    #Creating Labels for the plot
    aspectsLabel.append(''.join([l for l in aspect if l.isupper()]))
    
dfGenesEnrichments = pd.DataFrame(GOEnrichments, index=['Clusters of Genes','GO Annotations','Genes'], columns=aspectsLabel).transpose()
dfGenesEnrichments


# ### Visulize Enrichments

# In[ ]:


plotEnrichments(dfGenesEnrichments, dfDrugsEnrichments, saveFig)

