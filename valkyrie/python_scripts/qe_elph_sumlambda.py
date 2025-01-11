import os
import numpy as np

qpoint=16
qband=3
nsig=30
os.system('''grep Ef elph.inp_lambda.%s | awk '{print NR,$3,$8}' >DOS_Ef.txt''' % qpoint)
os.system('''grep lambda elph.inp_lambda.%s | awk '{print $3}' >tb.txt''' % qpoint)

qk_lambda=np.loadtxt('tb.txt')

sumlambda=[]
s=0
for ns in range(nsig):
    s=0
    for qb in range(qband):
        s+=qk_lambda[ns*qband+qb]
    sumlambda.append(s)
np.savetxt('lambda_vs_nsig'+str(qpoint)+'.dat',sumlambda)
print(sumlambda)
#os.system('rm tb.txt')

