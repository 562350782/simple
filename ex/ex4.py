import numpy as np

"---------------------------------------------Problem-------------------------------------------------------------------"

'两个同心圆柱形成的圆环，外圆柱半径为 0.1m，内圆柱体半径为0.01m' \
'圆环内有两种互不相融的不可压缩流体，两种流体的分界面为0.05m，流体密度分别为1000kg m和800kg m ，' \
'动力黏度分别为0.1Pas和0.01Pas。初始时刻，环内流体静止不动，t0时刻，两个同心圆柱分别以等角速度1rads和-1rads' \
'绕圆心柱o旋转。假设圆环内流体运动为层流流动。'

"------------------------------------------Initialization--------------------------------------------------------------"
r1 = 0.1   # 外径
r2 = 0.01  # 内径
r3 = 0.05  # 流体界面

n_in = 4 # 内流体划分
n_out = 5 # 外流体划分
omega1 = 1   # 外筒转动方向
omega2 = -1   # 内筒转动方向
detar_in = (r3-r2)/n_in  # 内流体的步长
detar_out = (r1-r3)/n_out  # 外流体的步长


rho1 = 1000.0  # 外流体密度
rho2 = 800.0  # 内流体密度
rho3 = (rho1 + rho2) / 2  # 两流体交界容积的密度

mu1 = 0.1  # 外流体动力粘度
mu2 = 0.01  # 外流体动力粘度
mu3 = (mu1 + mu2) / 2  # 两流体交界容积的粘度

A = np.zeros((n_in+n_out,n_in+n_out))  # n+1时层次的各项系数
u_old = np.zeros((n_in+n_out,1))  # n时层次的速度u
u_new = np.zeros((n_in+n_out,1))  # n+1时层次的速度u
r = np.zeros((n_in+n_out, 1))  # 创建数组保存每个节点的位置信息（r值）

"------------------------------------------Iterative_Solution----------------------------------------------------------"
'循环计算每个节点对应的位置信息（r值）'
for i in range(n_in+n_out):
    r[i] = (i+1/2)*detar_in + r2

'时间层走600次，每次0.1时长'
for t in range(600):
    detat = 0.1  # 设置时间层迭代的时间长度
    '高斯塞得尔迭代层，直接设置循环次数，经测试已经收敛'
    for j in range(20):
        '每个节点的计算层'
        for i in range (n_in+n_out): # python以0开始计算点，点的编号是0-89

            '内壁面边界条件设置主要参数有r2，omega2，rho2，mu2'
            if i == 0:
                'n+1时层速度u1的系数'
                A[i, i] = rho2 * detar_in / mu2 / detat + 2 / detar_in - 1 / (r[i + 1] + r[i]) + 1 / (r[i] + r2)
                'n+1时层速度u2的系数'
                A[i + 1, i] = 1 / detar_in + 1 / (r[i + 1] + r[i])
                'n+1时层内壁面边界处u的系数'
                boundary_in = 1 / detar_in - 1 / (r[i] + r2)
                'n时层内速度u1的系数'
                old = rho2 * detar_in / mu2 / detat
                'n+1时层速度u1的计算式'
                u_new[i] = (u_old[i, 0] * old + A[i + 1, i] * u_new[i + 1, 0] + boundary_in * omega2 * r2) / A[i, i]

            '内部流体节点速度的计算0-39'
            if i!= 0 and i < 4:
                'n+1时层速度ui的系数'
                A[i,i]=rho2 * detar_in/ mu2/detat + 2/detar_in - 1/(r[i+1]+r[i]) + 1/(r[i]+r[i-1])
                'n+1时层速度ui+1的系数'
                A[i+1, i] = 1/detar_in + 1/(r[i+1]+r[i])
                'n+1时层速度ui-1的系数'
                A[i-1, i] = 1/detar_in - 1/(r[i]+r[i-1])
                'n时层内速度u1的系数'
                old = rho2 * detar_in/ mu2/detat
                'n+1时层速度ui的计算式'
                u_new[i,0] =(u_old[i,0]*old + A[i+1, i]*u_new[i+1,0] + A[i-1, i]*u_new[i-1,0])/A[i,i]

            '内外流体边界处节点速度的计算'
            if i == 4:
                'n+1时层速度u40的系数'
                A[i, i] = rho3 * detar_in / mu3 / detat + 2 / detar_in - 1 / (r[i + 1] + r[i]) + 1 / (r[i] + r[i - 1])
                'n+1时层速度u41的系数'
                A[i + 1, i] = 1 / detar_in + 1 / (r[i + 1] + r[i])
                'n+1时层速度u39的系数'
                A[i - 1, i] = 1 / detar_in - 1 / (r[i] + r[i - 1])
                'n时层内速度u40的系数'
                old = rho3 * detar_in / mu3 / detat
                'n+1时层速度u40的计算式'
                u_new[i, 0] = (u_old[i, 0] * old + A[i + 1, i] * u_new[i + 1, 0] + A[i - 1, i] * u_new[i - 1, 0]) / A[i, i]

            '外部流体节点速度的计算41-88'
            if i > 4 and i!=8:
                'n+1时层速度ui的系数'
                A[i,i]=rho1 * detar_in / mu1 / detat + 2 / detar_in - 1/(r[i+1]+r[i]) + 1/(r[i]+r[i-1])
                'n+1时层速度ui+1的系数'
                A[i+1, i] = 1/detar_in + 1/(r[i+1]+r[i])
                'n+1时层速度ui-1的系数'
                A[i-1, i] = 1/detar_in - 1/(r[i]+r[i-1])
                'n时层内速度u1的系内壁面边界处u的系数数'
                old = rho1 * detar_in / mu1/detat
                'n+1时层速度ui的计算式'
                u_new[i,0] =(u_old[i , 0]*old + A[i+1, i]*u_new[i+1,0] + A[i-1, i]*u_new[i-1,0])/A[i,i]

            '外壁面边界条件设置主要参数有r1，omega1，rho1，mu1'
            if i == 8:
                'n+1时层速度u89的系数'
                A[i, i] = rho1 * detar_in / mu1 / detat +2 /detar_in +1 / (r1 + r[i]) + 1 / (r[i] + r[i - 1])
                'n+1时层内壁面边界处u的系数'
                boundary_out = 1/detar_in + 1 / (r1 + r[i])
                'n+1时层速度u88的系数'
                A[i - 1, i] = 1/detar_in - 1 / (r[i] + r[i - 1])
                'n时层内速度u89的系数'
                b = rho1 * detar_in / mu1 / detat
                'n+1时层速度u89的计算式'
                u_new[i] = (u_old[i ,0] * b + boundary_out * omega1 * r1 + A[i - 1, i] * u_new[i - 1, 0]) / A[i, i]

    '时层每走一层将速度保存用于下一时层的计算'
    u_old= u_new

"-------------------------------------------------result---------------------------------------------------------------"

'输出数据 格式：  节点位置  对应的速度'
for i in range(n_in+n_out):
    print('%f     %.15f' % (r[i] ,u_new[i]))