import numpy as np
from scipy import linalg
from panelvar.dynamic_panel_model import dynamic_panel_model

def PAR1_matrix(m: dynamic_panel_model):


    tbp = np.zeros((m.num_indep, m.num_indep), dtype='float64')
    if m.num_dep_lags==1:
        tbp[0:m.num_indep, 0:m.num_dep]=m.beta
    else:
        tbp[0:m.num_indep, 0:m.num_dep] = m.beta
        for i in range(0, m.num_indep-m.num_dep):
            tbp[i, i+m.num_dep]=1

    return tbp

def IRF_matrix(m: dynamic_panel_model, nahead: int):
    par1=PAR1_matrix(m)
    iden=np.identity(m.num_indep, dtype='float64')
    J=np.zeros((m.num_dep, m.num_indep),dtype='float64')
    for i in range(0,m.num_dep):
        J[i,i]=1
    # print('--J---')
    #
    #
    # print(J)
    list_tbr=[]
    list_tbr.append(np.identity(m.num_dep, dtype='float64'))
    temp=iden
    if nahead>=2:
        for i in range(1, nahead):
            temp=temp @ par1
            # print(temp)
            new_mat = J @ temp @ np.transpose(J)
            print("===new mat")
            print(new_mat)
            list_tbr.append(new_mat)

    return list_tbr

def oirf(m:dynamic_panel_model, nahead:int):
    ma_phi=IRF_matrix(m, nahead)
    print(m.residual.shape[0])
    print(m.num_obs)
    nona_residual=m.residual[~np.isnan(m.residual).any(axis=1)]
    Sigma_Hat1 = (np.cov(nona_residual, rowvar=False))
    print("sigma_hat1")
    print(Sigma_Hat1)
    p=linalg.cholesky(Sigma_Hat1)
    print("----chol")
    print(p)
    irf_output=[]

    MA_Phi_P=[]

    # Calculate 2.3.27 in Luetkepohl p. 58.
    for i0 in range(0,nahead):
        temp_mat=p @ ma_phi[i0]
        MA_Phi_P.append(temp_mat)
        #print("----------")
        #print(temp_mat)


    for i0 in range(0, m.num_dep):
        temp_mat=np.ndarray((nahead, p.shape[1]), dtype='float64')
        for i in range(0, nahead):
            temp_mat[i,:]=MA_Phi_P[i][i0,:]
        print(temp_mat)

    return irf_output
