#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


get_ipython().run_line_magic('run', '../scripts/ResultsAnalysis/enrichementAnalysis.py')


# ### VI

# In[ ]:


VI = pd.read_csv('../Data/GeneSets/VIsGenes.csv').values.flatten()
dfVIEnrichments = enrichmentAnalysisGprofiler(list(VI))


# In[ ]:


[t for t in dfVIEnrichments['term_name'].values if 'viral'in t]


# ### DEG

# In[ ]:


DEG = pd.read_csv('../Data/GeneSets/DEGsGenes.csv').values.flatten()
dfDEGEnrichments = enrichmentAnalysisGprofiler(list(DEG))


# In[ ]:


[t for t in dfDEGEnrichments['term_name'].values if 'viral' in t]


# ### VI-unique neighbors

# In[ ]:


VINeigh = pd.read_csv('../Data/GeneSets/VINeighboursUniqueGenes.csv').values.flatten()
dfVINeighEnrichments = enrichmentAnalysisGprofiler(list(VINeigh))


# In[ ]:


[t for t in dfVINeighEnrichments['term_name'].values if 'viral' in t]


# In[ ]:


dfVINeighEnrichments.sort_values(['source', 'p_value']).to_csv('../Data/Results/VINeighborsUnique_enrichment.csv', index=False)


# In[ ]:


VINeigh.shape, dfVINeighEnrichments.shape


# ### DEG-unique neighbors

# In[ ]:


DEGNeigh = pd.read_csv('../Data/GeneSets/DEGNeighboursUniqueGenes.csv').values.flatten()
dfDEGNeighnrichments = enrichmentAnalysisGprofiler(list(DEGNeigh))


# In[ ]:


[t for t in dfDEGNeighnrichments['term_name'].values if 'viral' in t]


# In[ ]:


dfDEGNeighnrichments.sort_values(['source', 'p_value']).to_csv('../Data/Results/DEGNeighborsUnique_enrichment.csv', index=False)


# In[ ]:


DEGNeigh.shape, dfDEGNeighnrichments.shape


# ### Background

# In[ ]:


dfTmp = pd.read_csv('../Data/GeneSets/allGenes.csv')
Background = dfTmp[dfTmp['Gene Description'] == 'Background']['Gene'].values
dfBackgroundEhnrichments = enrichmentAnalysisGprofiler(list(Background))


# In[ ]:


[t for t in dfBackgroundEhnrichments['term_name'].values if 'viral' in t]


# In[ ]:


dfBackgroundEhnrichments.sort_values(['source', 'p_value']).to_csv('../Data/Results/Background_enrichment.csv', index=False)


# In[ ]:


Background.shape, dfBackgroundEhnrichments.shape


# ### All Common Neighbors

# In[4]:


CN = pd.read_csv('../Data/GeneSets/CommonNeighboursGenes.csv').values.flatten()
dfCNEnrichments = enrichmentAnalysisGprofiler(list(CN))


# In[5]:


viralTerms = [t for t in dfCNEnrichments['term_name'].values if 'viral'in t]
dfCNEnrichments[dfCNEnrichments['term_name'].isin(viralTerms)][['p_value', 'source', 'term_name']]


# In[ ]:


dfCNEnrichments.to_csv('../Data/Results/CommonNeighbours_enrichment.csv', index=False)


# ##### Plot 

# In[20]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
dfCNEnrichmentsPlot = enrichmentAnalysisGprofiler(list(CN), toPlot=True)
dfCNEnrichmentsPlot.columns


# In[ ]:


get_ipython().run_cell_magic('R', '-i dfCNEnrichmentsPlot -w 8 -h 2 --units in -r 400', '# import df from global environment\n\n#Import libraries\nlibrary(enrichplot)\nlibrary(DOSE)\nlibrary(grid)\nlibrary(ggplot2)\n\n\n#Prepare the data for plotting\ngp_mod <- dfCNEnrichmentsPlot[,c("source", "term_id",\n                            "term_name", "p_value", "query_size",\n                            "intersection_size", "term_size",\n                            "effective_domain_size")]\ngp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)\ngp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)\nnames(gp_mod) = c("Category", "ID", "Description", "p.adjust",\n                  "query_size", "Count", "term_size", "effective_domain_size",\n                  "GeneRatio", "BgRatio")\nrow.names(gp_mod) = gp_mod$ID\n\n#Selecting the list of terms we want to plot\ntermsToPlot <- c(\'GO:0016032\', \'GO:0019083\', \'GO:0019080\', \'GO:0019058\', \'GO:0048524\')\n\n#Define as enrichResult object\ngp_mod_enrich = new("enrichResult", result = gp_mod[termsToPlot,])\n\n\n#Create plot\nbarplot(gp_mod_enrich, showCategory = 40, font.size = 8) + ggplot2::ylab("Intersection size")+ ggplot2::theme(plot.margin = unit(c(0,0,0,2), "cm"))\n\n#Save plot\nggsave("../Plots/Enrichments_CommonNeighbors_ViralProcesses.png")')


# In[ ]:





# In[ ]:





# ### Enrichment targeted common neighbours 

# #### Adding gene set information

# In[9]:


allGenes= pd.read_csv('../Data/GeneSets/allGenes.csv')
predictedDTIs = pd.read_csv('../Data/Results/predictedDTIs.csv')
dfDTIsGeneSet = pd.merge(predictedDTIs, allGenes, on='Gene')
dfDTIsGeneSet.columns = ['Gene', 'DrugBank ID', 'Score', 'Gene Description']
dfDTIsGeneSet


# #### Adding drugs information

# In[10]:


dfDrugGroup = pd.read_csv('../Data/Raw/Drug_Groups.csv')[['DrugBank ID', 'DrugBank Name', 'Group']]
dfDrugGroup.columns = ['DrugBank ID', 'DrugBank Name', 'Drug Status']
dfDTIsGeneSetDrugInfo = pd.merge(dfDTIsGeneSet, dfDrugGroup, on='DrugBank ID')

#Some drugs have different status, keeping the approved/experimental and if both are for the same drug keep approved
dfDTIsGeneSetDrugInfo = dfDTIsGeneSetDrugInfo[dfDTIsGeneSetDrugInfo['Drug Status'].isin(['approved', 'experimental'])]
indexDrugExperimentalAlsoApproved = dfDTIsGeneSetDrugInfo[dfDTIsGeneSetDrugInfo.duplicated(subset=['Gene', 'DrugBank ID'])].index
dfDTIsGeneSetDrugInfo.drop(indexDrugExperimentalAlsoApproved, axis=0, inplace=True)

dfDTIsGeneSetDrugInfo.to_csv('../Data/Results/predictedDTIsInfo.csv')


# In[11]:


dfDTIsCommonNeighApproved = dfDTIsGeneSetDrugInfo[(dfDTIsGeneSetDrugInfo['Drug Status']== 'approved')&(dfDTIsGeneSetDrugInfo['Gene Description']=='Common neighbor')]
dfDTIsCommonNeighApproved


# In[12]:


CNTargeted = dfDTIsCommonNeighApproved['Gene'].unique()
print(CNTargeted.shape[0], ' targeted common neighbors')
dfCNTargeted = enrichmentAnalysisGprofiler(list(CNTargeted))


# In[ ]:


dfCNTargeted = enrichmentAnalysisGprofiler(list(CNTargeted))
dfCNTargeted.to_csv('../Data/Results/TargetedCommonNeighbours_enrichment.csv', index=False)


# ##### Plot 

# In[ ]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
dfCNTargetedPlot = enrichmentAnalysisGprofiler(list(CNTargeted), toPlot=True)
dfCNTargetedPlot.columns


# In[ ]:


get_ipython().run_cell_magic('R', '-i dfCNTargetedPlot -w 8 -h 5 --units in -r 400', '# import df from global environment\n\n#Import libraries\nlibrary(enrichplot)\nlibrary(DOSE)\nlibrary(grid)\nlibrary(ggplot2)\n\n\n#Prepare the data for plotting\ngp_mod <- dfCNTargetedPlot[,c("source", "term_id",\n                            "term_name", "p_value", "query_size",\n                            "intersection_size", "term_size",\n                            "effective_domain_size")]\ngp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)\ngp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)\nnames(gp_mod) = c("Category", "ID", "Description", "p.adjust",\n                  "query_size", "Count", "term_size", "effective_domain_size",\n                  "GeneRatio", "BgRatio")\nrow.names(gp_mod) = gp_mod$ID\n\n#Selecting the list of terms we want to plot\ntermsToPlot <- c(\'GO:0007186\', \'GO:0098664\', \'GO:0007188\', \'GO:0007193\', \'GO:0007200\', \'GO:0008227\', \'GO:0004993\', \n                 \'GO:0000165\', \'GO:0043408\', \'GO:0043410\',    \n                 \'GO:0038086\', \'GO:0038091\',\n                 \'GO:0070371\', \'GO:0070372\', \'GO:0070374\', \n                 \'GO:0006198\', \'GO:0046058\', \'GO:0043951\', \'GO:0043949\', \'GO:0019933\', \'GO:0030552\',\n                 \'GO:0045834\', \'GO:0008015\', \'GO:0034702\', \'GO:0099589\', \'GO:0051378\', \'GO:0042417\', \'GO:1903351\',\n                 \'GO:1903350\', \'GO:0004517\',   \n                 \'KEGG:04151\', \'KEGG:04014\', \'KEGG:04010\', \'KEGG:04024\', \'KEGG:04370\',   \n                 \'REAC:R-HSA-375280\', \'REAC:R-HSA-1280215\', \'REAC:R-HSA-392154\')\n\n#Define as enrichResult object\ngp_mod_enrich = new("enrichResult", result = gp_mod[termsToPlot,])\n\n\n#Create plot\nbarplot(gp_mod_enrich, showCategory = 40, font.size = 8) + ggplot2::ylab("Intersection size")+ ggplot2::theme(plot.margin = unit(c(0,0,0,2), "cm"))\n\n#Save plot\nggsave("../Plots/Enrichments_targetedCommonNeighbors.png")')


# In[ ]:


dfCNTargeted.to_csv('../Data/Results/TargetedCommonNeighbours_enrichment.csv', index=False)

