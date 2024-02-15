#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd


# In[ ]:


get_ipython().run_line_magic('run', '../scripts/DataFusionFramework/choosingParameters.py')


# ### Creating matrix from data

# In[ ]:


proteinsNames = pd.read_csv('../Data/Names/ViralProteins.csv')['Proteins'].values
geneNames = pd.read_csv('../Data/Names/Genes.csv')['Genes'].values
drugNames = pd.read_csv('../Data/Names/Drugs.csv')['Drugs'].values
proteinsNames.shape[0], geneNames.shape[0], drugNames.shape[0]


# ### Execute computation of the dispersion coefficient for the different parametrizations

# In[ ]:


computeDispersionCoefficient([proteinsNames, geneNames, drugNames], path='../Data/Results/', savePath='../Data/Results/') 

