# import vaex as va
import pandas as pd
from  panelvar import regression



import time
start=time.time()
for i in range(1):

    df = pd.read_csv("data.csv")



    mydpd = regression.pvar(dep="n w k", df=df, lags=2, gmm="gmm(n w k, 2:.)",options="constant", identifiers=['id', 'year'], ahead=8, draw=200)
    print(mydpd.list_models[0].regression_result.beta)
    print(mydpd.list_models[0].regression_result.std_err)
#
#
    print(mydpd.list_models[0].irf)
print(time.time()-start)





