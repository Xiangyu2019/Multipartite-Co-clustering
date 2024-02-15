#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('run', '../scripts/Data/creatingMatrices.py')
get_ipython().run_line_magic('run', '../scripts/DataFusionFramework/framework.py')


# ### Creating matrix from data

# In[ ]:


r12Aug, l2, r23, l3, r12, a2, r23 = matricesCreation(saveMatrices=True, verbose=True, numpy=True, pathLoad='../Data/Raw/', pathSave='../Data/Matrices/', alpha=0.6)


# In[ ]:


r12Aug, l2, r23, l3 = matricesLoadNumpy(pathLoad='../Data/Matrices/', alpha=0.6)


# In[ ]:


r12Aug.shape, l2.shape, r23.shape, l3.shape


# In[ ]:


# Seeds used in the paper 
seeds = [3020, 3339, 2117, 4060, 2484, 4486, 1152, 3721,  656, 3287]


# ### Tunning parameters - compute the data fusion with different parameters

# In[ ]:


saveMatrices = True
verbose = True
loadPath='../Data/Matrices/'
savePath='../Data/Results/'
alpha = 0.6


for k1 in [3, 5]:
    for k2 in [80,100,60,120]:
        for k3 in [80,60,40]:

            #For each combination run 10 times the method
            for r, s in enumerate(seeds):

                if verbose: print('Parameters k1={}, k2={}, k3={}, alpha={}. Run number {}'.format(k1, k2, k3, alpha, r))
                G1, G2, G3, H12, H23 = matricesRandomAcolInitializationNumpy(r12Aug, r23, l3, k1, k2, k3, s)
                factors, errHistory = datafusionframework_triPMF(r12Aug, r23, l2, l3, G1, G2, G3, H12, H23,                                                           verbose=True, itStopChecking=100, epsilon=1e-5)
                H12, H23, G1, G2, G3 = factors

                #The folders must be created by hand (e.g. /0.6_3_80_100/)
                if saveMatrices:
                    np.save(savePath+'/{}_{}_{}_{}/G1_{}.npy'.format(str(alpha), str(k1), str(k2), str(k3), str(r)), G1)
                    np.save(savePath+'/{}_{}_{}_{}/G2_{}.npy'.format(str(alpha), str(k1), str(k2), str(k3), str(r)), G2)
                    np.save(savePath+'/{}_{}_{}_{}/G3_{}.npy'.format(str(alpha), str(k1), str(k2), str(k3), str(r)), G3)
                    np.save(savePath+'/{}_{}_{}_{}/H12_{}.npy'.format(str(alpha), str(k1), str(k2), str(k3), str(r)), H12)
                    np.save(savePath+'/{}_{}_{}_{}/H23_{}.npy'.format(str(alpha), str(k1), str(k2), str(k3), str(r)), H23)

                    np.savetxt(savePath+'/{}_{}_{}_{}/errHist_{}.npy'.format(str(alpha), str(k1), str(k2), str(k3), str(r)), errHistory)

