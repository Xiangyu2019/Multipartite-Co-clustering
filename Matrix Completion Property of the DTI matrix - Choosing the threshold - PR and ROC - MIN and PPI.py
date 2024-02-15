#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np


# In[ ]:


get_ipython().run_line_magic('run', '../scripts/ResultsAnalysis/PRandROCcurves.py')
get_ipython().run_line_magic('run', '../scripts/plots.py')


# ### MIN

# #### Load matrices

# In[ ]:


R23_reconstruct_MIN = reconstructR23('../Data/Results/')
r23_MIN = np.load('../Data/Matrices/Matrix_R23.npy')
r23_MIN.shape, R23_reconstruct_MIN.shape


# #### Find Threshold using PRcurve

# In[ ]:


maxF1_MIN, thresholdF1_MIN, precisionMaxF1_MIN, recallMaxF1_MIN, precision_MIN, recall_MIN, pr_auc_MIN = computePRCurve(r23_MIN.flatten(), R23_reconstruct_MIN.flatten())


# In[ ]:


saveFile = '../Plots/PR-curve.png'
plotPRCurve(precision_MIN, recall_MIN, maxF1_MIN, precisionMaxF1_MIN, recallMaxF1_MIN, pr_auc_MIN, thresholdF1_MIN, saveFile)


# #### ROC cruve

# In[ ]:


fpr_MIN, tpr_MIN, roc_auc_MIN = computeROC(r23_MIN.flatten(), R23_reconstruct_MIN.flatten())


# In[ ]:


saveFile = '../Plots/ROC-curve.png'
plotROCCurve(fpr_MIN, tpr_MIN,roc_auc_MIN, saveFile)


# #### Distribution of association scores

# In[ ]:


valuesNotInOriginal_MIN, alreadyPaired_inReconstructed_MIN = computeNewPairs(R23_reconstruct_MIN, r23_MIN)


# In[ ]:


saveFile = '../Plots/DTI_DistributionAssociationScores.png'
plotAssociationScores(valuesNotInOriginal_MIN, alreadyPaired_inReconstructed_MIN, thresholdF1_MIN, saveFile)


# #### Constrain the Reconstructed matrix --> predicted DTIs

# In[ ]:


geneNames = pd.read_csv('../Data/Names/Genes.csv')['Genes'].values
drugNames = pd.read_csv('../Data/Names/Drugs.csv')['Drugs'].values
geneNames.shape[0], drugNames.shape[0]


# In[ ]:


saveFile = '../Data/Results/predictedDTIs.csv'
restricDTIs(thresholdF1_MIN, R23_reconstruct_MIN, r23_MIN, drugNames, geneNames, saveFile)


# ## PPI

# #### Load Matrices

# In[ ]:


R23_reconstruct_PPI = reconstructR23('../Data/Results/PPIasHumanInteractome/')
r23_PPI = np.load('../Data/Matrices/PPIasHumanInteractome/Matrix_R23.npy')
r23_PPI.shape, R23_reconstruct_PPI.shape


# #### Find the threshold using PRcurve

# In[ ]:


maxF1_PPI, thresholdF1_PPI, precisionMaxF1_PPI, recallMaxF1_PPI, precision_PPI, recall_PPI, pr_auc_PPI = computePRCurve(r23_PPI.flatten(), R23_reconstruct_PPI.flatten())


# #### Plot PRcurve comparison MIN and PPI

# In[ ]:


saveFile = '../Plots/DTI_PRCurve_ComparisonMINPPI.png'
plotTwoPRCurve([precision_MIN, precision_PPI], 
                [recall_MIN, recall_PPI], 
                [maxF1_MIN, maxF1_PPI], 
                [thresholdF1_MIN, thresholdF1_PPI], 
                [precisionMaxF1_MIN, precisionMaxF1_PPI], 
                [recallMaxF1_MIN, recallMaxF1_PPI], 
                [pr_auc_MIN, pr_auc_PPI], saveFile)


# ### ROC

# In[ ]:


fpr_PPI, tpr_PPI, roc_auc_PPI = computeROC(r23_PPI.flatten(), R23_reconstruct_PPI.flatten())


# #### Plot ROCcurve comparison MIN and PPI

# In[ ]:


saveFile = '../Plots/DTI_ROCCurve_ComparisonMINPPI.png'
plotTwoROCcurves([fpr_MIN, fpr_PPI], [tpr_MIN, tpr_PPI], [roc_auc_MIN, roc_auc_PPI], saveFile)


# #### Distribution of association scores

# In[ ]:


valuesNotInOriginal_PPI, alreadyPaired_inReconstructed_PPI = computeNewPairs(R23_reconstruct_PPI, r23_PPI)


# In[ ]:


saveFile = '../Plots/DTI_DistributionAssociationScores_PPI.png'
plotAssociationScores(valuesNotInOriginal_PPI, alreadyPaired_inReconstructed_PPI, thresholdF1_PPI, saveFile)


# #### Constrain the Reconstructed matrix --> predicted DTIs

# In[ ]:


geneNames_PPI = pd.read_csv('../Data/Names/Genes_PPI.csv')['Genes'].values
drugNames = pd.read_csv('../Data/Names/Drugs.csv')['Drugs'].values
geneNames.shape[0], drugNames.shape[0]


# In[ ]:


saveFile = '../Data/Results/PPIasHumanInteractome/predictedDTIs_PPI.csv'
restricDTIs(thresholdF1_PPI, R23_reconstruct_PPI, r23_PPI, drugNames, geneNames_PPI, saveFile)

