from sklearn.model_selection import train_test_split
from sklearn import mixture
#from scikit_learn.sklearn import mixture
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import uuid
import os
from sklearn.decomposition import PCA
import random

nc = 100 #number of clusters
#run_id= str(uuid.uuid1())
run_id ="myoutput_lambda325_10b"
os.mkdir(run_id)

filename="testcontainer/mrm_sm_3500.csv"
driver_namesfile  ="testcontainer/drivers_promote.csv"
gene_namesfile = "testcontainer/targets_promote.csv"


#read files
driver_names =list(pd.read_csv(driver_namesfile, header=0, index_col=None).iloc[:,0])
gene_names = list(pd.read_csv(gene_namesfile, header=0, index_col=None).iloc[:, 0])
X_train_orig=pd.read_csv(filename, header=0, index_col=0)

#random seed
rsi=random.sample(range(1000),10)
#LAM=np.arange(250,500,50)
for i in range(1,11):

    idx=sorted(random.sample(range(46),int(round(46*0.8,0))))
    X_train=X_train_orig.iloc[:,idx]
    # PCA on data
    pca = PCA(n_components=X_train.values.shape[1])
    pca.fit(np.transpose(X_train.values))
    pca_score = pca.explained_variance_ratio_
    cumsum_pca_score=np.cumsum(pca_score)
    #keep 90% of variance
    num_components = np.where(cumsum_pca_score>0.9)[0][0]
    X_train.iloc[:, :]=np.transpose(pca.components_)

    #filter driver and gene matrix
    driver_matrix=X_train.loc[driver_names,:].iloc[: , 0:num_components]
    gene_matrix=X_train.loc[gene_names,:].iloc[:, 0:num_components]

    #center and scale
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    scaler.fit(np.transpose(gene_matrix))
    gene_matrix=np.transpose(scaler.transform(np.transpose(gene_matrix)))
    scaler.fit(np.transpose(driver_matrix))
    driver_matrix=np.transpose(scaler.transform(np.transpose(driver_matrix)))
    driver_matrix=driver_matrix.T

    print('Bootstrap '+ str(i))
    #print('lambda '+ str(LAM[i-1]))
    gmm = mixture.GaussianMixture(n_components=nc,
                              covariance_type="spherical", random_state=rsi[i-1], init_params="sparsGMM", max_iter=200, tol=10e-3, lambda_p=325)
    gmm.fit(gene_matrix, driver_matrix)
    labels = gmm.predict(gene_matrix, driver_matrix)
    prob_labels = gmm.predict_proba(gene_matrix, driver_matrix)
    labels_pd=pd.DataFrame(data=labels, index=gene_names, columns=np.arange(0, 1))
    prob_labels_pd = pd.DataFrame(data=prob_labels, index=gene_names, columns=np.arange(0, nc))
    betas_pd=pd.DataFrame(data=gmm.betas_,index=driver_names,columns=np.arange(0,nc))
    labels_pd.to_csv(run_id + "/labels" + "_" + str(i) +"_promote_"+ str(nc) +".csv", index=True, header=True)
    prob_labels_pd.to_csv(run_id + "/prob_labels" + "_" + str(i)+"_promote_"+ str(nc) + ".csv", index=True, header=True)
    betas_pd.to_csv(run_id + "/weights" + "_" + str(i) +"_promote_"+ str(nc)+ ".csv", index=True, header=True)
