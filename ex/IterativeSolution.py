import numpy as np

def solution(n,lamda,detax,T0,A,B,S):
    for i in range(n):
        if i == 0:
            A[i,i] = 3*lamda/detax
            A[i,i+1] = -lamda/detax
            B[i,0] = S*detax+2*lamda/detax*T0

        if i == n-1:
            A[i,i] = lamda/detax
            A[i,i-1] = -lamda/detax
            B[i,0] = S * detax

        if i != 0 and i != n-1:
            A[i,i] = 2*lamda/detax
            A[i,i-1] = -lamda/detax
            A[i,i+1] = -lamda/detax
            B[i,0]=S*detax

    x = np.linalg.solve(A,B)
    return x