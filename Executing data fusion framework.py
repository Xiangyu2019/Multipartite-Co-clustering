#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('run', '../scripts/Data/creatingMatrices.py')
get_ipython().run_line_magic('run', '../scripts/DataFusionFramework/framework.py')


# ### Creating/Loading matrices

# In[ ]:


PPIasHumanInteractome = True


# ##### Execute if the matrices need to be created

# In[ ]:


if PPIasHumanInteractome:
    r12Aug, l2, r23, l3, r12, a2, r23 = matricesCreationPPI(saveMatrices=True, verbose=True, numpy=True, pathLoad='../Data/Raw/', pathSave='../Data/Matrices/', alpha=0.6)
else: 
    r12Aug, l2, r23, l3, r12, a2, r23 = matricesCreation(saveMatrices=True, verbose=True, numpy=True, pathLoad='../Data/Raw/', pathSave='../Data/Matrices/', alpha=0.6)


# ##### Execute if the matrices are already created

# In[ ]:


if PPIasHumanInteractome:
    r12Aug, l2, r23, l3 = matricesLoadNumpyPPI(pathLoad='../Data/Matrices/', alpha=0.6)
else:
    r12Aug, l2, r23, l3 = matricesLoadNumpy(pathLoad='../Data/Matrices/', alpha=0.6)


# In[ ]:


r12Aug.shape, l2.shape, r23.shape, l3.shape


# ### Initializing matrices

# In[ ]:


#Parameters obtained after parametrization
k1, k2, k3 = 3, 120, 80


# In[ ]:


g1, g2, g3, h12, h23 = matricesSVDInitialization(r12Aug, r23, k1, k2, k3)


# In[ ]:


g1.shape, g2.shape, g3.shape, h12.shape, h23.shape 


# ### Execute data fusion

# In[ ]:


factors, errHistory = datafusionframework_triPMF(r12Aug, r23, l2, l3, g1, g2, g3, h12, h23, verbose=True, itStopChecking=100)

H12, H23, G1, G2, G3 = factors
np.save('../Data/Results/G1_SVD.npy', G1)
np.save('../Data/Results/G2_SVD.npy', G2)
np.save('../Data/Results/G3_SVD.npy', G3)
np.save('../Data/Results/H12_SVD.npy', H12)
np.save('../Data/Results/H23_SVD.npy', H23)

