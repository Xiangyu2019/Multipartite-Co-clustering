#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from scipy.stats import hypergeom, mannwhitneyu
import pandas as pd


# In[ ]:


get_ipython().run_line_magic('run', '../scripts/ResultsAnalysis/topologicalAnalysisNetwork.py')
get_ipython().run_line_magic('run', '../scripts/plots.py')


# ## Create Gene Sets

# #### Load VI and DEGs

# In[ ]:


VI = pd.read_csv('../Data/GeneSets/VIsGenes.csv').values.flatten()
DEG = pd.read_csv('../Data/GeneSets/DEGsGenes.csv').values.flatten()
print(VI.shape[0], ' VIs')
print(DEG.shape[0], 'DEGs')


# #### Compute gene sets for PPI

# In[ ]:


PPI = pd.read_csv('../Data/Matrices/PPIasHumanInteractome/Matrix_A2.csv', index_col=0)
computeGeneSets(PPI, VI, DEG, '../Data/GeneSets/PPIasHumanInteractome/')


# #### Compute gene sets for MIN

# In[ ]:


MIN = pd.read_csv('../Data/Matrices/Matrix_A2.csv', index_col=0)
computeGeneSets(MIN, VI, DEG, '../Data/GeneSets/')


# ## Overlap of the VI and DEGs neighbors 

# In[ ]:


DEGNeigh = pd.read_csv('../Data/GeneSets/DEGNeighboursGenes.csv').values.flatten()
VINeigh = pd.read_csv('../Data/GeneSets/VINeighboursGenes.csv').values.flatten()
DEGNeigh.shape[0], VINeigh.shape[0]


# In[ ]:


plotNeighborsOverlap(VINeigh, DEGNeigh)


# ###  Hypergeometric test
# 
# Let M be the total number of gene in the network, out of which K are in neigh(DEGs). N genes are in Neigh(VHIs), out of which X are also in neigh(DEGs). The probability of having X or more common neighbor by chance is obtained by the hypergeometric formula for enrichments.

# In[ ]:


commonNeigh = pd.read_csv('../Data/GeneSets/CommonNeighboursGenes.csv').values.flatten()
allGenes = pd.read_csv('../Data/Names/Genes.csv').values.flatten()


# In[ ]:


M = len(allGenes)
K = len(DEGNeigh)
N = len(VINeigh)
X = len(commonNeigh)
print('All genes ', M)
print('DEG neighbors ', K)
print('VI neighbros ', N)
print('Common neighbors ', X)
print('Percentage of common genes ', X/(K+N-X)*100)
print('p-value: ', hypergeom.sf(X-1, M, K, N))


# ### Degree 

# In[ ]:


MIN = pd.read_csv('../Data/Matrices/Matrix_A2.csv', index_col=0)


# In[ ]:


allGenes = pd.read_csv('../Data/GeneSets/allGenes.csv')
allGenes.columns = ['Gene', 'Gene Set']


# In[ ]:


degrees = computeDegree(MIN)
dfDegreeGeneSet = pd.merge(degrees, allGenes, on='Gene')


# In[ ]:


saveFile = 'Plot/BoxPlotComparisonDegrees.png'
plotBoxPlotDegreeComparison(dfDegreeGeneSet, saveFile=None)


# In[ ]:


dfCentralityMeasures = centralityMeasures(MIN, verbose=True)
dfCentralityMeasures.to_csv('../Data/Results/CentralityMeasuresTmp.csv')


# In[ ]:


dfCentralityMeasuresGeneSet = pd.merge(dfCentralityMeasures, allGenes, on='Gene')
dfCentralityMeasuresGeneSet.groupby('Gene Set').mean().to_csv('../Data/Results/CentralityMeasuresMeans.csv')


# #### Test for statistically significant 
# 
#  Two-sided Mann-Whitney-Wilcoxon test (p < 0.05).

# In[ ]:


res = []
geneSets = dfDegreeGeneSet['Gene Set'].unique()
for i, id1 in enumerate(geneSets):
    geneset1 = dfDegreeGeneSet[dfDegreeGeneSet['Gene Set'] == id1]['Degree'].values
    for id2 in list(geneSets)[i+1:]:
        geneset2 = dfDegreeGeneSet[dfDegreeGeneSet['Gene Set'] == id2]['Degree'].values
        manmannwhitneyuRes = mannwhitneyu(geneset1, geneset2, use_continuity=False, alternative='two-sided')
        res.append([id1, id2, manmannwhitneyuRes[0], manmannwhitneyuRes[1]])
                    
                    
dfRes = pd.DataFrame(res, columns=['Gene Set 1','Gene Set 2', 'u-statistic', 'p-value'])                    
dfRes.to_csv('../Data/Results/MolecularNetworkAnalysis-DegreeComparison.csv', index=False)


# ## GDV

# ### Compare GDV of the different gene sets in the MIN

# #### Load gene sets

# In[ ]:


dfAllGenes = pd.read_csv('../Data/GeneSets/allGenes.csv')
VI = dfAllGenes[dfAllGenes['Gene Description'] == 'VI']['Gene'].unique()
DEG = dfAllGenes[dfAllGenes['Gene Description'] == 'DEG']['Gene'].unique()
VINeigh = dfAllGenes[dfAllGenes['Gene Description'] == 'VI-unique neighbor']['Gene'].unique()
DEGNeigh = dfAllGenes[dfAllGenes['Gene Description'] == 'DEG-unique neighbor']['Gene'].unique()
ComNeigh = dfAllGenes[dfAllGenes['Gene Description'] == 'Common neighbor']['Gene'].unique()
Background = dfAllGenes[dfAllGenes['Gene Description'] == 'Background']['Gene'].unique()

print(VI.shape[0], ' VIs')
print(DEG.shape[0], ' DEGs')
print(VINeigh.shape[0], ' VI-unique Neighbors')
print(DEGNeigh.shape[0], ' DEG-unique Neighbors')
print(ComNeigh.shape[0], ' Common Neighbors')
print(Background.shape[0], ' Background')


# #### Load MIN

# In[ ]:


MIN = pd.read_csv('../Data/Matrices/Matrix_A2.csv', index_col=0)


# #### Create orca files

# In[ ]:


fileName = "MIN_names_edgelist"
label_mapping, reverse_mapping = createOrcaFile(MIN, fileName)


# #### Plot

# In[ ]:


geneSets = [VI, DEG, VINeigh, DEGNeigh, ComNeigh, Background]
nameSubGeneSets = ['VI', 'DEG', 'VI-unique neighbors', 'DEG-unique neighbors', 'Common neighbors', 'Background']
colors = ["#7fb9d7", "#fdb456", "#1f78b4", "#ff7f00", '#33a02c', "#a4a4a4"]
saveFile = '../Plots/GDV_signature_geneGroups.jpg'
outFilePath = '../Data/OrcaFormatFiles'
plotGDVDifferentSetsSameNetwork(fileName, outFilePath, reverse_mapping, geneSets, nameSubGeneSets, colors, saveFile)

