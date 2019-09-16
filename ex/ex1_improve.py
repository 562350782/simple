import IterativeSolution
import Initialization
import numpy as np
"长度为1m，横截面积为0.01m2的圆柱型金属棒，一端横截面露出，其余部分被厚的绝热橡胶层包裹。金属棒表面均匀缠绕着热电阻丝，热电阻丝的发热总功率为20W"
"将横截面露出的一端放入100℃沸腾的水池中，已知该金属材料的导热系数平均为 5w/mk，放置很长时间后求解金属棒的温度分布"
"------------------------------------------Initialization--------------------------------------------------------------"
lamda= 5
S=2000
n=50
l=1
detax=l/n
print(detax)
A = np.zeros((50, 50))
B = np.zeros((50, 1))
T = np.zeros((50, 1))
"------------------------------------------Iterative_Solution----------------------------------------------------------"
# def solution(x,A,B,lamda,detax,S):
for j in range(50000):
    for i in range(n):
        if i == 0:
            A[i, i + 1] = lamda / detax
            A[i, i] = 3 * lamda / detax
            B[i, 0] = S*detax + 200 * lamda / detax
            T[i, 0] = (B[i, 0] + A[i, i + 1]*T[i + 1, 0])/A[i, i]
            # print(A)
        if i == n-1:
            A[i, i] = lamda / detax
            A[i, i - 1] = lamda / detax
            B[i, 0] = S*detax
            T[i, 0] = (B[i, 0] + A[i, i - 1]*T[i - 1, 0])/A[i, i]
            # print(T[i, 0])
        if i != 0 and i != n - 1:
            A[i, i] = 2 * lamda/detax
            A[i, i - 1] = lamda/detax
            A[i, i + 1] = lamda/detax
            B[i, 0] = S*detax
            T[i, 0] = (B[i, 0] + A[i, i-1]*T[i-1, 0] + A[i, i+1]*T[i+1, 0])/A[i, i]
print(T)
"------------------------------------------main--------------------------------------------------------------"
# result = solution(x,A,B,lamda,detax,S)
# slove = IterativeSolution.solution(Initialization.n, Initialization.lamda, Initialization.detax,
#                                        Initialization.T0, Initialization.A, Initialization.B, Initialization.S)
# print(slove)

