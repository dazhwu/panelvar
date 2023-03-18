# import vaex as va
import pandas as pd
from  pypvar import regression



import time
start=time.time()
for i in range(1):

    df = pd.read_csv("data.csv")


    #command_str='y L(1:?).y L(1:?).x  | gmm(y, 2:3) iv(L(1:1).x)| timedumm'
    #command_str='n L1.n |gmm(n, 2:4)'
    #mydpd = regression.abond(command_str, df, ['id', 'year'])
    # df = pd.read_csv("data.csv")
    #command_str='n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k)  '
    #command_str = 'n L(1:?).n w L1.w L2.w L(3:?).w k  | gmm(n, 2:3) pred(w k) |hqic  timedumm'
    #command_str="n  L1.n  L2.n  w  k | gmm(n, 2:3) pred(w k) "
    #command_str="n  L1.n  w  L1.w  L2.w  L3.w  L4.w  k | gmm(n, 2:3) pred(w k) |hqic"

    #pvar(dep, df: DataFrame, identifiers: list, exog = [], lags = 1, gmm = "", options_str = "" )
   # mydpd = regression.pvar(dep="n w k", df=df, lags=0, gmm="gmm(n w k, 2:?)",options="collapse", identifiers=['id', 'year'])
    mydpd = regression.pvar(dep="n w", df=df, lags=2, gmm="gmm(n w k, 2:.)",options="", identifiers=['id', 'year'], ahead=8, draw=200)
# print(mydpd.list_models[0].regression_result.beta)
# print(mydpd.list_models[0].regression_result.std_err)
#
#
# print(mydpd.list_models[0].irf)
print(time.time()-start)





