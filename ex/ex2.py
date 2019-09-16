import numpy as np

"------------------------------------------Initialization--------------------------------------------------------------"
T_up = 2
T_down = 1
Tf = 1
hf = 10
q = 1
S = 1
lamda = 1
detax = 1
A = np.zeros((3,3))
T = np.zeros((3,3))
# print(A)
# print(T)
"------------------------------------------Iterative_Solution----------------------------------------------------------"
for k in range(5000):
    for j in range(3):
        for i in range(3):
            "===================================第一行y = 0（j = 0）的迭代========================================="
            "当P为(0,0)点时"
            if j == 0 and i == 0:
                A[i, j] = 22 / 3 * lamda / detax
                A[i+1, j] = 1 * lamda / detax
                A[i, j+1] = 1 * lamda / detax
                B = 19 / 3 * S
                T[i, j]  = (B + A[i+1, j]*T[i+1, j] + A[i, j+1]*T[i, j+1])/A[i, j]
            if j == 0 and i != 0 and i != 2:
                A[i, j] = 5 * lamda / detax
                A[i + 1, j] = 1 * lamda / detax
                A[i, j + 1] = 1 * lamda / detax
                A[i - 1, j] = 1 * lamda / detax
                B = 3 * S
                T[i, j] = (B + A[i + 1, j] * T[i + 1, j] + A[i, j + 1] * T[i, j + 1] + A[i - 1, j] * T[i - 1, j]) / A[i, j]
            if j == 0 and i == 2:
                A[i, j] = 4 * lamda / detax
                A[i, j + 1] = 1 * lamda / detax
                A[i - 1, j] = 1 * lamda / detax
                B = 5 * S
                T[i, j] = (B + A[i, j + 1] * T[i, j + 1] + A[i - 1, j] * T[i - 1, j]) / A[i, j]
            "===================================第一行y = 1（j = 1）的迭代========================================="
            if j == 1 and i == 0:
                A[i, j] = 19 / 3 * lamda / detax
                A[i + 1, j] = 1 * lamda / detax
                A[i, j + 1] = 1 * lamda / detax
                A[i, j - 1] = 1 * lamda / detax
                B = 13 / 3 * S
                T[i, j] = (B + A[i + 1, j] * T[i + 1, j] + A[i, j + 1] * T[i, j + 1] + A[i, j - 1] * T[i, j - 1]) / A[i, j]
            if j == 1 and i != 0 and i != 2:
                A[i, j] = 4 * lamda / detax
                A[i + 1, j] = 1 * lamda / detax
                A[i, j + 1] = 1 * lamda / detax
                A[i - 1, j] = 1 * lamda / detax
                A[i, j - 1] = 1 * lamda / detax
                B = 1 * S
                T[i, j] = (B + A[i + 1, j] * T[i + 1, j] + A[i, j + 1] * T[i, j + 1] + A[i - 1, j] * T[i - 1, j] + A[i, j - 1] * T[i, j - 1]) / A[i, j]
            if j == 1 and i == 2:
                A[i, j] = 3 * lamda / detax
                A[i, j + 1] = 1 * lamda / detax
                A[i - 1, j] = 1 * lamda / detax
                A[i, j - 1] = 1 * lamda / detax
                B = 3 * S
                T[i, j] = (B + A[i, j + 1] * T[i, j + 1] + A[i - 1, j] * T[i - 1, j] + A[i, j - 1] * T[i, j - 1]) / A[i, j]
            "===================================第一行y = 2（j = 2）的迭代========================================="
            if j == 2 and i == 0 :
                A[i, j] = 22 / 3 * lamda / detax
                A[i + 1, j] = 1 * lamda / detax
                A[i, j - 1] = 1 * lamda / detax
                B = 25 / 3 * S
                T[i, j] = (B + A[i + 1, j] * T[i + 1, j] + A[i, j - 1] * T[i, j - 1]) / A[i, j]
            if j == 2 and i != 0 and i != 2:
                A[i, j] = 5 * lamda / detax
                A[i + 1, j] = 1 * lamda / detax
                A[i - 1, j] = 1 * lamda / detax
                A[i, j - 1] = 1 * lamda / detax
                B = 5 * S
                T[i, j] = (B + A[i + 1, j] * T[i + 1, j] + A[i - 1, j] * T[i - 1, j] + A[i, j - 1] * T[i, j - 1]) / A[i, j]
            if j == 2 and i == 2:
                A[i, j] = 4 * lamda / detax
                A[i - 1, j] = 1 * lamda / detax
                A[i, j - 1] = 1 * lamda / detax
                B = 7 * S
                T[i, j] = (B + A[i - 1, j] * T[i - 1, j] + A[i, j - 1] * T[i, j - 1]) / A[i, j]

"-------------------------------------------------result---------------------------------------------------------------"
# print(A)
print(T)