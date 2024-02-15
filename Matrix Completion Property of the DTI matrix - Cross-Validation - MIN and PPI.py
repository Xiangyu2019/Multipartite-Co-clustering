#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_score, recall_score, f1_score, precision_recall_fscore_support
import pandas as pd
import numpy as np
import random


# In[ ]:


get_ipython().run_line_magic('run', '../scripts/Data/creatingMatrices.py')
get_ipython().run_line_magic('run', '../scripts/DataFusionFramework/framework.py')
get_ipython().run_line_magic('run', '../scripts/ResultsAnalysis/PRandROCcurves.py')
get_ipython().run_line_magic('run', '../scripts/plots.py')


# In[ ]:


PPIasHumanInteractome = False


# In[ ]:


if PPIasHumanInteractome:
    pathResults = '../Data/Results/PPIasHumanInteractome/'
    pathPlots = '../Plots/PPIasHumanInteractome/'
    pathMatrices = '../Data/Matrices/PPIasHumanInteractome/'
else:
    pathResults = '../Data/Results/'
    pathPlots = '../Plots/'
    pathMatrices = '../Data/Matrices/'


# ### Loading matrices

# In[ ]:


if PPIasHumanInteractome:
    r12Aug, l2, r23, l3 = matricesLoadNumpyPPI(pathLoad='../Data/Matrices/', alpha=0.6)
else:
    r12Aug, l2, r23, l3 = matricesLoadNumpy(pathLoad='../Data/Matrices/', alpha=0.6)


# In[ ]:


r12Aug.shape, l2.shape, r23.shape, l3.shape


# In[ ]:


knownDTIs = np.where(r23 == 1)
knownDTIsList = np.concatenate((knownDTIs[0], knownDTIs[1])).reshape((2,knownDTIs[0].shape[0])).transpose()

notKnownDTIs = np.where(r23 == 0)
notKnownDTIsList = np.concatenate((notKnownDTIs[0], notKnownDTIs[1])).reshape((2,notKnownDTIs[0].shape[0])).transpose()

DTIsList = np.concatenate((knownDTIsList, notKnownDTIsList))
DTIsLabels = np.concatenate((np.ones(knownDTIsList.shape[0]), np.zeros(notKnownDTIsList.shape[0])))

del knownDTIs
del knownDTIsList
del notKnownDTIs
del notKnownDTIsList


# ### Cross-validation r23 (DTIs) in 10-folds

# In[ ]:


verbose=True

pr_aucs_train, pr_aucs_validation = [], []
roc_aucs_train, roc_aucs_validation = [], []

precisions_train, precisions_validation = [], []
recalls_train, recalls_validation = [], []

fprs_train, fprs_validation = [], []
tprs_train, tprs_validation = [], []

skf = StratifiedKFold(n_splits=10, random_state=42, shuffle=True)
fold = 0
for train_index, validation_index in skf.split(DTIsList, DTIsLabels):
    fold += 1

    if verbose: print("FOLD", fold)
    if verbose: print("  TRAIN:", len(train_index))
    if verbose: print("  TEST:", len(validation_index))

    trainingknownDTIs = DTIsList[train_index][DTIsLabels[train_index] == 1]

    #Construct DTIs matrix with training entries
    r23_training_DTIs = np.zeros(r23.shape)
    r23_training_DTIs[trainingknownDTIs[:, 0], trainingknownDTIs[:, 1]] = 1
    if verbose: print(r23_training_DTIs.shape, r23_training_DTIs.sum())


    #Initializing factor matrices
    k1, k2, k3 = 3, 120, 80 #Parameters obtained after parametrization
    g1, g2, g3, h12, h23 = matricesSVDInitialization(r12Aug, r23_training_DTIs, k1, k2, k3)
    if verbose: print('  Factors shapes', g1.shape, g2.shape, g3.shape, h12.shape, h23.shape)

        
    #Computing decoposition  
    if verbose: print('  ->Staring decompostion')
    factors, errHistory = datafusionframework_triPMF(r12Aug, r23_training_DTIs, l2, l3, g1, g2, g3, h12, h23, verbose=verbose, itStopChecking=100)
    del errHistory
    
    H12, H23, G1, G2, G3 = factors
    np.save(pathResults + 'CrossFoldValidation/DecompositionFactors/G1_SVD_10Fold_{}.npy'.format(fold), G1)
    np.save(pathResults + 'CrossFoldValidation/DecompositionFactors/G2_SVD_10Fold_{}.npy'.format(fold), G2)
    np.save(pathResults + 'CrossFoldValidation/DecompositionFactors/G3_SVD_10Fold_{}.npy'.format(fold), G3)
    np.save(pathResults + 'CrossFoldValidation/DecompositionFactors/H12_SVD_10Fold_{}.npy'.format(fold), H12)
    np.save(pathResults + 'CrossFoldValidation/DecompositionFactors/H23_SVD_10Fold_{}.npy'.format(fold), H23)

    #Relese memory
    del G1, H12
    del g1, g2, g3, h12, h23 
    
    #Reconstruct matrix
    R23_reconstruct = np.matmul(G2, np.matmul(H23, np.transpose(G3)))
    #Relese memory
    del G2, G3, H23, factors
    
    #Separe reconstructed matrix in train and validation
    R23_reconstruct_training = R23_reconstruct[DTIsList[train_index][:,0], DTIsList[train_index][:,1]]
    R23_reconstruct_validating = R23_reconstruct[DTIsList[validation_index][:,0], DTIsList[validation_index][:,1]]
    #Relese memory
    del R23_reconstruct

    #Separe input matrix in train and validation
    r23_training = r23[DTIsList[train_index][:,0], DTIsList[train_index][:,1]]
    r23_validating = r23[DTIsList[validation_index][:,0], DTIsList[validation_index][:,1]]


    #Computing PR-Curves
    if verbose: print('  ->Computing PR-Curves')
    precision_train, recall_train, pr_auc_train = computePRCurve(r23_training, R23_reconstruct_training, computeThreshold=False)
    precision_validation, recall_validation, pr_auc_validation = computePRCurve(r23_validating, R23_reconstruct_validating, computeThreshold=False)
    
    pr_aucs_train.append(pr_auc_train)
    pr_aucs_validation.append(pr_auc_validation)
    
    precisions_train.append(precision_train)
    precisions_validation.append(precision_validation)
    
    recalls_train.append(recall_train)
    recalls_validation.append(recall_validation)
    

    #Computing ROC-Curves
    if verbose: print('  ->Computing ROC-Curves')
    fpr_train, tpr_train, roc_auc_train = computeROC(r23_training, R23_reconstruct_training)
    fpr_validation, tpr_validation, roc_auc_validation = computeROC(r23_validating, R23_reconstruct_validating)    

    roc_aucs_train.append(roc_auc_train)
    roc_aucs_validation.append(roc_auc_validation)
    
    fprs_train.append(fpr_train)
    fprs_validation.append(fpr_validation)
    
    tprs_train.append(tpr_train)
    tprs_validation.append(tpr_validation)

    
del R23_reconstruct_training, R23_reconstruct_validating
    
if verbose: print('Saving Results')
#Saving results
fileSaveAUCs = pathResults + 'CrossFoldValidation/AUCsResults.csv'
dfCVResults = pd.DataFrame([pr_aucs_train, pr_aucs_validation, roc_aucs_train, roc_aucs_validation], 
                           columns=['Fold 1', 'Fold 2', 'Fold 3', 'Fold 4', 'Fold 5', 'Fold 6', 'Fold 7', 'Fold 8', 'Fold 9', 'Fold 10'],
                           index=['PR_AUC_Train', 'PR_AUC_Validation', 'ROC_AUC_Train', 'ROC_AUC_Validation']).to_csv(fileSaveAUCs)

np.save(pathResults + 'CrossFoldValidation/precision-train.npy', precisions_train)
np.save(pathResults + 'CrossFoldValidation/precision-validation.npy', precisions_validation)
np.save(pathResults + 'CrossFoldValidation/recall-train.npy', recalls_train)
np.save(pathResults + 'CrossFoldValidation/recall-validation.npy', recalls_validation)
np.save(pathResults + 'CrossFoldValidation/thresholds-train.npy', threshold_train)
np.save(pathResults + 'CrossFoldValidation/thresholds-validation.npy', threshold_validation)
np.save(pathResults + 'CrossFoldValidation/fpr-train.npy', fprs_train)
np.save(pathResults + 'CrossFoldValidation/fpr-validation.npy', fprs_validation)
np.save(pathResults + 'CrossFoldValidation/tpr-train.npy', tprs_train)
np.save(pathResults + 'CrossFoldValidation/tpr-validation.npy', tprs_validation)


# ### Plot PR and ROC curves

# In[ ]:


if 'pr_aucs_train' not in locals():
    dfCVResults = pd.read_csv(pathResults + 'CrossFoldValidation/AUCsResults.csv', index_col=0)
    pr_aucs_train = dfCVResults.loc['PR_AUC_Train'].values
    pr_aucs_validation = dfCVResults.loc['PR_AUC_Validation'].values
    roc_aucs_train = dfCVResults.loc['ROC_AUC_Train'].values
    roc_aucs_validation = dfCVResults.loc['ROC_AUC_Validation'].values


# In[ ]:


if 'precisions_train' not in locals():
    precisions_train = np.load(pathResults + 'CrossFoldValidation/precision-train.npy', allow_pickle=True)
if 'recalls_train' not in locals():
    recalls_train = np.load(pathResults + 'CrossFoldValidation/recall-train.npy', allow_pickle=True)
plotCVPRCurve(precisions_train, recalls_train, pr_aucs_train, titleExtra=' (Train)', saveFile=pathPlots + 'CrossFoldValidation/PR-Curve-Train.png')
del precisions_train, recalls_train


# In[ ]:


if 'precisions_validation' not in locals():    
    precisions_validation = np.load(pathResults + 'CrossFoldValidation/precision-validation.npy', allow_pickle=True)
if 'recalls_validation' not in locals():
    recalls_validation = np.load(pathResults + 'CrossFoldValidation/recall-validation.npy', allow_pickle=True)
plotCVPRCurve(precisions_validation, recalls_validation, pr_aucs_validation, titleExtra=' (Validation)', saveFile='../Plots/CrossFoldValidation/PR-Curve-Validation.png')
del precisions_validation, recalls_validation


# In[ ]:


if 'fprs_train' not in locals():
    fprs_train = np.load(pathResults + 'CrossFoldValidation/fpr-train.npy', allow_pickle=True)
if 'tprs_train' not in locals():
    tprs_train = np.load(pathResults + 'CrossFoldValidation/tpr-train.npy', allow_pickle=True)
plotCVROCcurves(fprs_train, tprs_train, roc_aucs_train,  titleExtra=' (Train)', saveFile=pathPlots + 'CrossFoldValidation/ROC-Curve-Train.png')
del fprs_train, tprs_train


# In[ ]:


if 'fprs_validation' not in locals():
    fprs_validation = np.load(pathResults + 'CrossFoldValidation/fpr-validation.npy', allow_pickle=True)
if 'tprs_validation' not in locals():
    tprs_validation = np.load(pathResults + 'CrossFoldValidation/tpr-validation.npy', allow_pickle=True)
plotCVROCcurves(fprs_validation, tprs_validation, roc_aucs_validation,  titleExtra=' (Validation)', saveFile= pathPlots + 'CrossFoldValidation/ROC-Curve-Validation.png')
del fprs_validation, tprs_validation


# ## Compute precision and recall for the choosen threshold

# In[ ]:


verbose = False
t = 0.295

precisionsThreshold, recallsThreshold = [], []

skf = StratifiedKFold(n_splits=10, random_state=42, shuffle=True)
fold = 0
for train_index, validation_index in skf.split(DTIsList, DTIsLabels):
    fold += 1
    print("FOLD", fold)

    trainingknownDTIs = DTIsList[train_index][DTIsLabels[train_index] == 1]

    #Construct DTIs matrix with training entries
    r23_training_DTIs = np.zeros(r23.shape)
    r23_training_DTIs[trainingknownDTIs[:, 0], trainingknownDTIs[:, 1]] = 1
    if verbose: print(r23_training_DTIs.shape, r23_training_DTIs.sum())

    G2 = np.load(pathResults + 'CrossFoldValidation/DecompositionFactors/G2_SVD_10Fold_{}.npy'.format(fold))
    G3 = np.load(pathResults + 'CrossFoldValidation/DecompositionFactors/G3_SVD_10Fold_{}.npy'.format(fold))
    H23 = np.load(pathResults + 'CrossFoldValidation/DecompositionFactors/H23_SVD_10Fold_{}.npy'.format(fold))

    
    #Reconstruct matrix
    R23_reconstruct = np.matmul(G2, np.matmul(H23, np.transpose(G3)))
    
    #Separe reconstructed matrix in train and validation
    R23_reconstruct_training = R23_reconstruct[DTIsList[train_index][:,0], DTIsList[train_index][:,1]]
    R23_reconstruct_validating = R23_reconstruct[DTIsList[validation_index][:,0], DTIsList[validation_index][:,1]]
    #Relese memory
    del R23_reconstruct
    
    R23_reconstruct_training_threshold = np.where(R23_reconstruct_training > t, 1, 0)
    R23_reconstruct_validating_threshold = np.where(R23_reconstruct_validating > t, 1, 0)
    
    #Separate input matrix in train and validation
    r23_training = r23[DTIsList[train_index][:,0], DTIsList[train_index][:,1]]
    r23_validating = r23[DTIsList[validation_index][:,0], DTIsList[validation_index][:,1]]
    
    
    #Compute Precision and Recall
    precisionsThreshold.append(precision_score(r23_validating,R23_reconstruct_validating_threshold))
    recallsThreshold.append(recall_score(r23_validating,R23_reconstruct_validating_threshold))
    
np.save(pathResults + 'CrossFoldValidation/precisionThreshold-validation.npy', precisionsThreshold)
np.save(pathResults + 'CrossFoldValidation/recallThreshold-validation.npy', recallsThreshold)


# In[ ]:


if 'pr_aucs_train' not in locals():
    dfCVResults = pd.read_csv(pathResults + 'CrossFoldValidation/AUCsResults.csv', index_col=0)
    pr_aucs_validation = dfCVResults.loc['PR_AUC_Validation'].values

if 'precisions_validation' not in locals():    
    precisions_validation = np.load(pathResults + 'CrossFoldValidation/precision-validation.npy', allow_pickle=True)
if 'recalls_validation' not in locals():
    recalls_validation = np.load(pathResults + 'CrossFoldValidation/recall-validation.npy', allow_pickle=True)
if 'precisionsThreshold' not in locals():
    precisionsThreshold = np.load(pathResults + 'CrossFoldValidation/precisionThreshold-validation.npy')
if 'recallsThreshold' not in locals():    
    recallsThreshold = np.load(pathResults + 'CrossFoldValidation/recallThreshold-validation.npy')

    
plotCVPRCurve(precisions_validation, recalls_validation, pr_aucs_validation, thresholdPerformance=[precisionsThreshold, recallsThreshold], titleExtra=' (Validation)', saveFile= pathPlots + 'CrossFoldValidation/PR-Curve-Validation-Threshold.png')
del pr_aucs_validation, dfCVResults, precisions_validation, recalls_validation, precisionsThreshold, recallsThreshold


# In[ ]:




