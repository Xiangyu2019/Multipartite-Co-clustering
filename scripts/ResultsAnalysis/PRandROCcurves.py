from sklearn.metrics import precision_recall_curve, average_precision_score, roc_auc_score, roc_curve, roc_curve, auc

import pandas as pd
import numpy as np


def reconstructR23(savePath):
    """
    Reconstructs the R23 with the matrix factors from the decomposition
    
    Parameters
    ----------
    
    savePath : string
        Indicates the path where the matrix factors are stored
        
    Return
    ------
    Returns a numpy array with th values of the recontructed R23

    """
    
    G2 = np.load(savePath+'G2_SVD.npy')
    G3 = np.load(savePath+'G3_SVD.npy')
    H23 = np.load(savePath+'H23_SVD.npy')
    
    R23_reconstruct = np.matmul(G2, np.matmul(H23, np.transpose(G3)))

    return R23_reconstruct


def computePRCurve(gT, predicted, computeThreshold=True):
    """
    Computes the precision, recall and fbeta-score for each threshold.  
    
    Parameters
    ----------
    
    gT : numpy array
        Contains the ground truth values (e.g. original R23)
    predicted : numpy array
        Contains the  predicted values (e.g., reconstructed R23)
    computeThreshold: boolean
        Indicates if the threshold corresponding to the maximum in fbeta-score must be computed.
    Return
    ------
    Returns the maximum fbeta-score, the threshold, prescision and recall associated to the maxim fbeta-score; the precision and recall values for each threhsold and the arear under the curve.

    """
    

    precision, recall, thresholdsPR = precision_recall_curve(gT, predicted)
    pr_auc = auc(recall,precision)

    if computeThreshold:
        f1 = [2 * (p * r) / (p + r) for p, r in zip(precision,recall)]

        pMaxF1 = np.argmax(np.array(f1))   

        return max(f1), thresholdsPR[pMaxF1], precision[pMaxF1], recall[pMaxF1], precision, recall, pr_auc
    else:
        return precision, recall, pr_auc

def computeROC(gT, predicted):
    """
    Computes the fpr and tprs for each threshold
    
    Parameters
    ----------
    
    gT : numpy array
        Contains the ground truth values (e.g. original R23)
    predicted : numpy array
        Contains the  predicted values (e.g., reconstructed R23)
        
    Return
    ------
    Returns the false positive rate and true positive rate for each threshold; and the arear under the curve.

    """
            
    fpr, tpr, thresholdsROC = roc_curve(gT, predicted)
    roc_auc = roc_auc_score(gT, predicted)

        
    return fpr, tpr, roc_auc


def computeNewPairs(R23_reconstruct, r23, verbose=True):
    """
    Computes the scores of those pairs not in the original DTI and the scores of those values in the original DTI
    
    Parameters
    ----------
    
    r23 : numpy array
        Contains the original R23
    R23_reconstruct : numpy array
        Contains the reconstructed R23
    verbose: boolean
        Indicates whether information of the number of pairs have to be printed. By default True
        
    Return
    ------
    Returns the scores of those pairs in the original DTI and not in the original DTI separatelly.

    """
    
    valuesNotInOriginal = R23_reconstruct[np.where(r23 == 1.0, False, True)]

    #Obtain original paired values
    row_org, col_org = np.where(r23 == 1.0)
    org = [[r, c] for r, c in zip(row_org, col_org)]
    if verbose: print(len(org), 'original Paired')

    #Obtain values of the alreadyPairedValues
    alreadyPaired_inReconstructed = [R23_reconstruct[op_row, op_col] for op_row, op_col in org]
    if verbose: print(len(alreadyPaired_inReconstructed), 'already paired in the recontructed')
        
    return valuesNotInOriginal, alreadyPaired_inReconstructed


def restricDTIs(thresholdF1, R23_reconstruct, r23, drugNames, geneNames, saveFile, verbose=True):
    """
    Computes the scores of those pairs not in the original DTI and the scores of those values in the original DTI
    
    Parameters
    ----------
    thresholdF1 : float
        Indicates the threshold used to restric the DTIs
    R23_reconstruct : numpy array
        Contains the reconstructed R23
    r23 : numpy array
        Contains the original R23
    drugNames : list
        Contains the name of the drugs
    geneNames : list
        Contain the name of the genes
    saveFile : string
        Indicates the path and name of the file where to save the predicted DTIs list.
    verbose: boolean
        Indicates whether information of the new DTIs have to be printed. By default True
        
    Return
    ------
    Returns the scores of those pairs in the original DTI and not in the original DTI separatelly.

    """
        
    DF_R23_reconstruct = pd.DataFrame(R23_reconstruct)
    df = DF_R23_reconstruct[DF_R23_reconstruct > thresholdF1]

    row_org, col_org = np.where(r23 == 1.0)
    org = [[r, c] for r, c in zip(row_org, col_org)]

    #Remove the orginal pairs
    dfNotOrg = df.copy()
    count = 0
    for row, col in org:
        if not np.isnan(dfNotOrg.iloc[row, col]): count += 1
        dfNotOrg.iloc[row, col] = np.nan


    #Create the list of new pairs
    dfNotOrg.columns = drugNames
    dfNotOrg.index = geneNames
    dfDTIs = dfNotOrg.reset_index().melt(id_vars=['index'], value_vars=dfNotOrg.columns).dropna().sort_values('value', ascending=False)
    dfDTIs.columns = ['Gene', 'Drug', 'Score']
    dfDTIs.reset_index(inplace=True,drop=True)
    if verbose: print('With threshold {} ({}), the number of new pairs is {}; {} of original pairs'.format(thresholdF1, 'F1', dfDTIs.shape[0], count))
    dfDTIs.to_csv(saveFile, index=False)
