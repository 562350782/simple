import IterativeSolution
import Initialization
"长度为1m，横截面积为0.01m2的圆柱型金属棒，一端横截面露出，其余部分被厚的绝热橡胶层包裹。金属棒表面均匀缠绕着热电阻丝，热电阻丝的发热总功率为20W"
"将横截面露出的一端放入100℃沸腾的水池中，已知该金属材料的导热系数平均为 5w/mk，放置很长时间后求解金属棒的温度分布"
"------------------------------------------Initialization--------------------------------------------------------------"
# lamda= 5
# S=2000
# n=50
# l=1.0
# detax=l/n
# print(detax)
# A = np.zeros((50,50))
# B=np.zeros((50,1))
# T0 = 100
"------------------------------------------Iterative_Solution----------------------------------------------------------"
# def solution(n,lamda,detax,T0):
#     for i in range(n):
#         if i == 0:
#             A[i,i] = 3*lamda/detax
#             A[i,i+1] = -lamda/detax
#             B[i,0] = S*detax+2*lamda/detax*T0
#             print(A)
#             print(B)
#         if i == n-1:
#             A[i,i] = lamda/detax
#             A[i,i-1] = -lamda/detax
#             B[i,0] = S * detax
#             print(A)
#             print(B)
#         if i != 0 and i != n-1:
#             A[i,i] = 2*lamda/detax
#             A[i,i-1] = -lamda/detax
#             A[i,i+1] = -lamda/detax
#             B[i,0]=S*detax
#             print(A)
#             print(B)
#     x = np.linalg.solve(A,B)
#     return x
"------------------------------------------main--------------------------------------------------------------"
def main():
    slove = IterativeSolution.solution(Initialization.n, Initialization.lamda, Initialization.detax,
                                       Initialization.T0, Initialization.A, Initialization.B, Initialization.S)
    print(slove)

if __name__ == '__main__':
    main()

