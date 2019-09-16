import numpy as np
"---------------------------------------------Problem-------------------------------------------------------------------"

'假设0＜=x,y<=1的方腔内充满粘性不可压缩流体，左右，下壁固定，上壁以u=1运动，试求Re=100,200,400时的定常解'

"------------------------------------------Initialization--------------------------------------------------------------"
l = 1  # x,y的长度
n = 15  # 内节点法,x,y方向各有１０个节点
detar_x = l / n
detar_y = l / n
detar_t = 1
Re = 100
m = 0

show_u = np.zeros(((n+1)*(n+1), 3))
show_v = np.zeros(((n+1)*(n+1), 3))
k = 0

a_u = np.zeros((n+1, n+1))  # n+1时层u的各项系数
a_uu = np.zeros((n+1, n+1))  # n时层u的各项系数
a_u_add = np.zeros((n+1, n+1))

u = np.zeros((n+1, n+1))  # n+1时层次的速度u
uu = np.zeros((n+1, n+1))  # n时层次的速度u
u_add = np.zeros((n+1, n+1))

v = np.zeros((n+1, n+1))  # n时层次的速度u
vv = np.zeros((n+1, n+1))  # n+1时层次的速度u
v_add = np.zeros((n+1, n+1))

a_p = np.zeros((n+1, n+1))
a_e = np.zeros((n+1, n+1))
a_w = np.zeros((n+1, n+1))
a_n = np.zeros((n+1, n+1))
a_s = np.zeros((n+1, n+1))

p = np.zeros((n+1, n+1))
pp = np.zeros((n+1, n+1))
p_add = np.zeros((n+1, n+1))

d = np.zeros((n+1, n+1))

i = 1
j = 1
i_u = 1
j_u = 1
i_v = 1
j_v = 1
E = 1
W = 1
N = 1
S = 1

'亚松弛'
alpha_u = 0.7
alpha_v = 0.7
alpha_p = 0.7
"------------------------------------------Iterative_Solution----------------------------------------------------------"
u[n] = 1


def u_solution(x, y, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u):
    '设置detar_x与detar_y的长度，主要用于边界点和角点，当为左右边界１．５detar_x，当为上下边界1.5detar_y'
    '左下角点'
    if x == 1 and y == 0:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

        a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
        a_n[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
        a_w[y, x] = 0
        a_s[y, x] = 0

    '左上角点'
    if x == 1 and y == n - 1:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

        a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
        a_n[y, x] = 0
        a_w[y, x] = 0
        a_s[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

    '右下角点'
    if x == n - 1 and y == 0:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

        a_e[y, x] = 0
        a_n[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
        a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
        a_s[y, x] = 0

    '右上角点'
    if x == n - 1 and y == n - 1:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

        a_e[y, x] = 0
        a_n[y, x] = 0
        a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
        a_s[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

    '上边界点'
    if y == n - 1 and x != 0 and x != 1 and x != n - 1 and x != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * detar_x / detar_y

        a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
        a_n[y, x] = 0
        a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
        a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

    '下边界点'
    if y == 0 and x != 0 and x != 1 and x != n - 1 and x != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * detar_x / detar_y

        a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
        a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
        a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
        a_s[y, x] = 0

    '右边界点'
    if x == n - 1 and y != n - 1 and y != 0 and y != n:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x - 1]+ uu[y, x])/2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * 1.5 * detar_x / detar_y

        a_e[y, x] = 0
        a_n[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
        a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
        a_s[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

    '左边界点'
    if x == 1 and y != n - 1 and y != 0 and y != n:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * 1.5 * detar_x / detar_y

        a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
        a_n[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
        a_w[y, x] = 0
        a_s[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

    '内点'
    if x != 0 and x != n - 1 and x != n and y != 0 and y != n - 1 and y != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
        a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
        a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
        a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

    b = -detar_x * detar_y / detar_t * uu[y, x]
    P = detar_y * (p[y, x] - p[y, x - 1])

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    # a_e[y, x] = a_e[y, x] * E
    # a_n[y, x] = a_n[y, x] * N
    # a_w[y, x] = a_w[y, x] * W
    # a_s[y, x] = a_s[y, x] * S

    if x == 1 and y == n - 1:
        b += -2 * 1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
    if x == n - 1 and y == n - 1:
        b += -2 * 1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
    if y == n - 1 and x != 0 and x != 1 and x != n - 1 and x != n:
        b += -2 * detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)

    '亚松弛设置,注意两者的顺序'
    b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
    a_p[y, x] = a_p[y, x] / alpha_u

    u[y, x] = (a_e[y, x] * u[y, x + 1]
               + a_n[y, x] * u[y + 1, x]
               + a_w[y, x] * u[y, x - 1]
               + a_s[y, x] * u[y - 1, x]
               + b + P) / (-a_p[y, x])

    return u[y, x]

def v_solution(x, y, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_v):
    '设置detar_x与detar_y的长度，主要用于边界点和角点，当为左右边界１．５detar_x，当为上下边界1.5detar_y'
    '左下角点'
    if x == 0 and y == 1:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-vv[y - 1, x], 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0)
        a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
        a_w[y, x] = 0
        a_s[y, x] = 0

    '左上角点'
    if x == 0 and y == n - 1:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0)
        a_n[y, x] = 0
        a_w[y, x] = 0
        a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

    '右下角点'
    if x == n - 1 and y == 1:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = 0
        a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
        a_w[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
        a_s[y, x] = 0

    '右上角点'
    if x == n - 1 and y == n - 1:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = 0
        a_n[y, x] = 0
        a_w[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
        a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

    '左边界'
    if x == 0 and y != 0 and y != 1 and y != n - 1 and y != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0)
        a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
        a_w[y, x] = 0
        a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

    '右边界'
    if x == n - 1 and y != 0 and y != 1 and y != n - 1 and y != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = 0
        a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
        a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
        a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

    '下边界'
    if y == 1 and x != 0 and x != n - 1 and x != n:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-vv[y - 1, x], 0)) \
                    + 2 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0)
        a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
        a_w[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
        a_s[y, x] = 0

    '上边界'
    if y == n - 1 and x != 0 and x != n - 1 and x != n:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 2 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0)
        a_n[y, x] = 0
        a_w[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
        a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

    '内节点'
    if x != 0 and x != n - 1 and x != n and y != 0 and y != 1 and y != n - 1 and y != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0)
        a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
        a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[y-1, x] + uu[y, x]) / 2, 0)
        a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

    b = -detar_x * detar_y / detar_t * vv[y, x]
    P = detar_x * (p[y, x] - p[y - 1, x])

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    # a_e[y, x] = a_e[y, x] * E
    # a_n[y, x] = a_n[y, x] * N
    # a_w[y, x] = a_w[y, x] * W
    # a_s[y, x] = a_s[y, x] * S

    '亚松弛设置'
    b = b + (1 - alpha_v) * a_p[y, x] / alpha_v * vv[y, x]
    a_p[y, x] = a_p[y, x] / alpha_v

    v[y, x] = (a_e[y, x] * v[y, x + 1]
               + a_n[y, x] * v[y + 1, x]
               + a_w[y, x] * v[y, x - 1]
               + a_s[y, x] * v[y - 1, x]
               + b + P) / (-a_p[y, x])

    return v[y, x]

def p_add_solution(x, y, i_u, j_u, i_v, j_v, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p):
    # detar_x_u = i_u * detar_x
    # detar_y_u = j_u * detar_y
    #
    # detar_x_v = i_v * detar_x
    # detar_y_v = j_v * detar_y

    '左上角点'
    if x==0 and y == n-1:
        a_e[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

        a_w[y, x] = 0

        a_n[y, x] = 0

        a_s[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y
    '左下角点'
    if x == 0 and y == 0:
        a_e[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

        a_w[y, x] = 0

        a_n[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_s[y, x] = 0
    '右上角点'
    if x == n-1 and y == n-1:
        a_e[y, x] = 0

        a_w[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

        a_n[y, x] = 0

        a_s[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '右下角点'
    if x == n-1 and y == 0:
        a_e[y, x] = 0

        a_w[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

        a_n[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_s[y, x] = 0

    '上边界点'
    if y == n-1 and x != 0 and x != n-1:
        a_e[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * detar_x / detar_y

        a_w[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * detar_x / detar_y

        a_n[y, x] = 0

        a_s[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 2 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '下边界点'
    if y == 0 and x != 0 and x != n - 1:
        a_e[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * detar_x / detar_y

        a_w[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * detar_x / detar_y

        a_n[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                    + 2 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_s[y, x] = 0

    '左边界点'
    if x == 0 and y != 0 and y != n-1:
        a_e[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * 1.5 * detar_x / detar_y

        a_w[y, x] = 0

        a_n[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                    + 3 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_s[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '右边界点'
    if x == n-1 and y != 0 and y != n - 1:
        a_e[y, x] = 0

        a_w[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * 1.5 * detar_x / detar_y

        a_n[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                    + 3 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_s[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '内节点'
    if y != 0 and y != n-1 and x != 0 and x != n-1:
        a_e[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x+1] + vv[y, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_w[y, x] = detar_x * detar_y / detar_t \
                        + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                        + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                        + 2 / Re * detar_y / detar_x \
                        + 2 / Re * detar_x / detar_y

        a_n[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y , x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

        a_s[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    d[y, x + 1] = detar_y * alpha_u / a_e[y, x]
    d[y, x - 1] = detar_y * alpha_u / a_w[y, x]
    d[y + 1, x] = detar_x * alpha_v / a_n[y, x]
    d[y - 1, x] = detar_x * alpha_v / a_s[y, x]

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    # d[y, x + 1] = d[y, x + 1] * E
    # d[y + 1, x] = d[y + 1, x] * N
    # d[y, x - 1] = d[y, x - 1] * W
    # d[y - 1, x] = d[y - 1, x] * S

    a_p[y, x + 1] = d[y, x + 1] * detar_y
    a_p[y, x - 1] = d[y, x - 1] * detar_y
    a_p[y + 1, x] = d[y + 1, x] * detar_x
    a_p[y - 1, x] = d[y - 1, x] * detar_x

    a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
    b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

    p_add[y, x] = (
                      a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1]
                      + a_p[y + 1, x] * p_add[y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add
                  ) / a_p[y, x]

    p[y, x] = pp[y, x] + alpha_p * p_add[y, x]

    return p[y, x]

def u_add_solution (x, y, i_u, j_u, i_v, j_v, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p):
    detar_x_u = i_u * detar_x
    detar_y_u = j_u * detar_y

    '左下角点'
    if x == 1 and y == 0:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

    '左上角点'
    if x == 1 and y == n - 1:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

    '右下角点'
    if x == n - 1 and y == 0:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

    '右上角点'
    if x == n - 1 and y == n - 1:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * 1.5 * detar_x / detar_y

    '上边界点'
    if y == n - 1 and x != 0 and x != 1 and x != n - 1 and x != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * detar_x / detar_y

    '下边界点'
    if y == 0 and x != 0 and x != 1 and x != n - 1 and x != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 3 / Re * detar_x / detar_y

    '右边界点'
    if x == n - 1 and y != n - 1 and y != 0 and y != n:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x - 1]+ uu[y, x])/2, 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * 1.5 * detar_x / detar_y

    '左边界点'
    if x == 1 and y != n - 1 and y != 0 and y != n:
        a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                    + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * 1.5 * detar_x / detar_y

    '内点'
    if x != 0 and x != n - 1 and x != n and y != 0 and y != n - 1 and y != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    d[y, x - 1] = detar_y * alpha_u / a_p[y, x]

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    # d[y, x - 1] = d[y, x - 1] * W

    u[y, x] = u[y, x] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])

    return u[y, x]

def v_add_solution (x, y, i_u, j_u, i_v, j_v, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p):
    detar_x_v = i_v * detar_x
    detar_y_v = j_v * detar_y

    '左下角点'
    if x == 0 and y == 1:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-vv[y - 1, x], 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '左上角点'
    if x == 0 and y == n - 1:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '右下角点'
    if x == n - 1 and y == 1:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '右上角点'
    if x == n - 1 and y == n - 1:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '左边界'
    if x == 0 and y != 0 and y != 1 and y != n - 1 and y != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '右边界'
    if x == n - 1 and y != 0 and y != 1 and y != n - 1 and y != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 3 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '下边界'
    if y == 1 and x != 0 and x != n - 1 and x != n:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-vv[y - 1, x], 0)) \
                    + 2 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '上边界'
    if y == n - 1 and x != 0 and x != n - 1 and x != n:
        a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                    + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 2 / Re * 1.5 * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    '内节点'
    if x != 0 and x != n - 1 and x != n and y != 0 and y != 1 and y != n - 1 and y != n:
        a_p[y, x] = detar_x * detar_y / detar_t \
                    + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                    + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                    + 2 / Re * detar_y / detar_x \
                    + 2 / Re * detar_x / detar_y

    d[y - 1, x] = detar_x * alpha_v / a_p[y, x]

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    # d[y - 1, x] = d[y - 1, x] * S

    v[y, x] = v[y, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

    return v[y, x]

for k in range(100):
    "------------------------------------------u控制方程计算---------------------------------------------"
    for y in range(n):
        for x in range(1, n):
            '左下角点'
            if x == 1 and y == 0:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

                a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
                a_n[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2,
                                                                                0)
                a_w[y, x] = 0
                a_s[y, x] = 0

                b = -1.5 * detar_x * detar_y / detar_t * uu[y, x]
                P = detar_y * 1.5 * (p[y, x] - p[y, x - 1])

            '左上角点'
            if x == 1 and y == n - 1:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

                a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
                a_n[y, x] = 0
                a_w[y, x] = 0
                a_s[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

                b = -1.5 * detar_x * detar_y / detar_t * uu[y, x]
                P = detar_y * 1.5 * (p[y, x] - p[y, x - 1])

            '右下角点'
            if x == n - 1 and y == 0:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

                a_e[y, x] = 0
                a_n[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2,
                                                                                0)
                a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
                a_s[y, x] = 0

                b = -1.5 * detar_x * detar_y / detar_t * uu[y, x]
                P = detar_y * 1.5 * (p[y, x] - p[y, x - 1])

            '右上角点'
            if x == n - 1 and y == n - 1:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

                a_e[y, x] = 0
                a_n[y, x] = 0
                a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
                a_s[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

                b = -1.5 * detar_x * detar_y / detar_t * uu[y, x]
                P = detar_y * 1.5 * (p[y, x] - p[y, x - 1])

            '上边界点'
            if y == n - 1 and x != 0 and x != 1 and x != n - 1 and x != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * detar_x / detar_y

                a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
                a_n[y, x] = 0
                a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
                a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

                b = -detar_x * detar_y / detar_t * uu[y, x]
                P = detar_y * (p[y, x] - p[y, x - 1])

            '下边界点'
            if y == 0 and x != 0 and x != 1 and x != n - 1 and x != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * detar_x / detar_y

                a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
                a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
                a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
                a_s[y, x] = 0

                b = -detar_x * detar_y / detar_t * uu[y, x]
                P = detar_y * (p[y, x] - p[y, x - 1])

            '右边界点'
            if x == n - 1 and y != n - 1 and y != 0 and y != n:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x - 1] + uu[y, x]) / 2, 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 1.5 * detar_x / detar_y

                a_e[y, x] = 0
                a_n[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2,
                                                                                0)
                a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
                a_s[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

                b = -1.5 * detar_x * detar_y / detar_t * uu[y, x]
                P = detar_y * 1.5 * (p[y, x] - p[y, x - 1])

            '左边界点'
            if x == 1 and y != n - 1 and y != 0 and y != n:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 1.5 * detar_x / detar_y

                a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
                a_n[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2,
                                                                                0)
                a_w[y, x] = 0
                a_s[y, x] = -1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

                b = -1.5 * detar_x * detar_y / detar_t * uu[y, x]
                P = detar_y * 1.5 * (p[y, x] - p[y, x - 1])

            '内点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != n - 1 and y != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)
                a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
                a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[x - 1, y] + uu[y, x]) / 2, 0)
                a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x - 1] + vv[y, x]) / 2, 0)

                b = -detar_x * detar_y / detar_t * uu[y, x]
                P = detar_y * (p[y, x] - p[y, x - 1])

            if x == 1 and y == n - 1:
                b += -2 * 1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
            if x == n - 1 and y == n - 1:
                b += -2 * 1.5 * detar_x / detar_y / Re - 1.5 * detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)
            if y == n - 1 and x != 0 and x != 1 and x != n - 1 and x != n:
                b += -2 * detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0)

            '亚松弛设置,注意两者的顺序'
            b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
            a_p[y, x] = a_p[y, x] / alpha_u

            u[y, x] = (a_e[y, x] * u[y, x + 1]
                       + a_n[y, x] * u[y + 1, x]
                       + a_w[y, x] * u[y, x - 1]
                       + a_s[y, x] * u[y - 1, x]
                       + b + P) / (-a_p[y, x])

    "------------------------------------------v控制方程计算---------------------------------------------"
    for y in range(1, n):
        for x in range(n):
            '左下角点'
            if x == 0 and y == 1:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-vv[y - 1, x], 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2,
                                                                                0)
                a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
                a_w[y, x] = 0
                a_s[y, x] = 0

                b = -detar_x * 1.5 * detar_y / detar_t * vv[y, x]
                P = detar_x * 1.5 * (p[y, x] - p[y - 1, x])

            '左上角点'
            if x == 0 and y == n - 1:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2,
                                                                                0)
                a_n[y, x] = 0
                a_w[y, x] = 0
                a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

                b = -detar_x * 1.5 * detar_y / detar_t * vv[y, x]
                P = detar_x * 1.5 * (p[y, x] - p[y - 1, x])

            '右下角点'
            if x == n - 1 and y == 1:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = 0
                a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
                a_w[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
                a_s[y, x] = 0

                b = -detar_x * 1.5 * detar_y / detar_t * vv[y, x]
                P = detar_x * 1.5 * (p[y, x] - p[y - 1, x])

            '右上角点'
            if x == n - 1 and y == n - 1:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = 0
                a_n[y, x] = 0
                a_w[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
                a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

                b = -detar_x * 1.5 * detar_y / detar_t * vv[y, x]
                P = detar_x * 1.5 * (p[y, x] - p[y - 1, x])

            '左边界'
            if x == 0 and y != 0 and y != 1 and y != n - 1 and y != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0)
                a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
                a_w[y, x] = 0
                a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

                b = -detar_x * detar_y / detar_t * vv[y, x]
                P = detar_x * (p[y, x] - p[y - 1, x])
            '右边界'
            if x == n - 1 and y != 0 and y != 1 and y != n - 1 and y != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = 0
                a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
                a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
                a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

                b = -detar_x * detar_y / detar_t * vv[y, x]
                P = detar_x * (p[y, x] - p[y - 1, x])

            '下边界'
            if y == 1 and x != 0 and x != n - 1 and x != n:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-vv[y - 1, x], 0)) \
                            + 2 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2,
                                                                                0)
                a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
                a_w[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
                a_s[y, x] = 0

                b = -detar_x * 1.5 * detar_y / detar_t * vv[y, x]
                P = detar_x * 1.5 * (p[y, x] - p[y - 1, x])

            '上边界'
            if y == n - 1 and x != 0 and x != n - 1 and x != n:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 2 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2,
                                                                                0)
                a_n[y, x] = 0
                a_w[y, x] = -1.5 * detar_y / detar_x / Re - 1.5 * detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
                a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

                b = -detar_x * 1.5 * detar_y / detar_t * vv[y, x]
                P = detar_x * 1.5 * (p[y, x] - p[y - 1, x])

            '内节点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != 1 and y != n - 1 and y != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-(uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0)
                a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)
                a_w[y, x] = -detar_y / detar_x / Re - detar_y * max((uu[y - 1, x] + uu[y, x]) / 2, 0)
                a_s[y, x] = -detar_x / detar_y / Re - detar_x * max((vv[y, x] + vv[y - 1, x]) / 2, 0)

                b = -detar_x * detar_y / detar_t * vv[y, x]
                P = detar_x * (p[y, x] - p[y - 1, x])

            '亚松弛设置'
            b = b + (1 - alpha_v) * a_p[y, x] / alpha_v * vv[y, x]
            a_p[y, x] = a_p[y, x] / alpha_v

            v[y, x] = (a_e[y, x] * v[y, x + 1]
                       + a_n[y, x] * v[y + 1, x]
                       + a_w[y, x] * v[y, x - 1]
                       + a_s[y, x] * v[y - 1, x]
                       + b + P) / (-a_p[y, x])

    "------------------------------------------压力修正---------------------------------------------"
    for y in range(n):
        for x in range(n):
            '左上角点'
            if x == 0 and y == n - 1:
                a_e[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                            + 1.5 * detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

                a_w[y, x] = 0

                a_n[y, x] = 0

                a_s[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                d[y, x + 1] = detar_y * alpha_u / a_e[y, x]
                d[y, x - 1] = 0
                d[y + 1, x] = 0
                d[y - 1, x] = detar_x * alpha_v / a_s[y, x]

            '左下角点'
            if x == 0 and y == 0:
                a_e[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                            + 1.5 * detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

                a_w[y, x] = 0

                a_n[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_s[y, x] = 0

                d[y, x + 1] = detar_y * alpha_u / a_e[y, x]
                d[y, x - 1] = 0
                d[y + 1, x] = detar_x * alpha_v / a_n[y, x]
                d[y - 1, x] = 0

            '右上角点'
            if x == n - 1 and y == n - 1:
                a_e[y, x] = 0

                a_w[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

                a_n[y, x] = 0

                a_s[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                d[y, x + 1] = 0
                d[y, x - 1] = detar_y * alpha_u / a_w[y, x]
                d[y + 1, x] = 0
                d[y - 1, x] = detar_x * alpha_v / a_s[y, x]

            '右下角点'
            if x == n - 1 and y == 0:
                a_e[y, x] = 0

                a_w[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

                a_n[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_s[y, x] = 0

                d[y, x + 1] = 0
                d[y, x - 1] = detar_y * alpha_u / a_w[y, x]
                d[y + 1, x] = detar_x * alpha_v / a_n[y, x]
                d[y - 1, x] = 0

            '上边界点'
            if y == n - 1 and x != 0 and x != n - 1:
                a_e[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * detar_x / detar_y

                a_w[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * detar_x / detar_y

                a_n[y, x] = 0

                a_s[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 2 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                d[y, x + 1] = detar_y * alpha_u / a_e[y, x]
                d[y, x - 1] = detar_y * alpha_u / a_w[y, x]
                d[y + 1, x] = 0
                d[y - 1, x] = detar_x * alpha_v / a_s[y, x]

            '下边界点'
            if y == 0 and x != 0 and x != n - 1:
                a_e[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * detar_x / detar_y

                a_w[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * detar_x / detar_y

                a_n[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                            + 2 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_s[y, x] = 0

                d[y, x + 1] = detar_y * alpha_u / a_e[y, x]
                d[y, x - 1] = detar_y * alpha_u / a_w[y, x]
                d[y + 1, x] = detar_x * alpha_v / a_n[y, x]
                d[y - 1, x] = 0

            '左边界点'
            if x == 0 and y != 0 and y != n - 1:
                a_e[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                            + 1.5 * detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 1.5 * detar_x / detar_y

                a_w[y, x] = 0

                a_n[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                            + 3 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_s[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                d[y, x + 1] = detar_y * alpha_u / a_e[y, x]
                d[y, x - 1] = 0
                d[y + 1, x] = detar_x * alpha_v / a_n[y, x]
                d[y - 1, x] = detar_x * alpha_v / a_s[y, x]

            '右边界点'
            if x == n - 1 and y != 0 and y != n - 1:
                a_e[y, x] = 0

                a_w[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + 1.5 * detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 1.5 * detar_x / detar_y

                a_n[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                            + 3 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_s[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                d[y, x + 1] = 0
                d[y, x - 1] = detar_y * alpha_u / a_w[y, x]
                d[y + 1, x] = detar_x * alpha_v / a_n[y, x]
                d[y - 1, x] = detar_x * alpha_v / a_s[y, x]

            '内节点'
            if y != 0 and y != n - 1 and x != 0 and x != n - 1:
                a_e[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x + 1] + uu[y, x + 2]) / 2, 0) + max(-(uu[y, x + 1] + uu[y, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x + 1] + vv[y + 1, x]) / 2, 0) + max(-(vv[y, x + 1] + vv[y, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_w[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_n[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x + 1] + uu[y + 1, x + 1]) / 2, 0) + max(-(uu[y + 1, x] + uu[y, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 2, x] + vv[y + 1, x]) / 2, 0) + max(-(vv[y + 1, x] + vv[y, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                a_s[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

                d[y, x + 1] = detar_y * alpha_u / a_e[y, x]
                d[y, x - 1] = detar_y * alpha_u / a_w[y, x]
                d[y + 1, x] = detar_x * alpha_v / a_n[y, x]
                d[y - 1, x] = detar_x * alpha_v / a_s[y, x]

            a_p[y, x + 1] = d[y, x + 1] * detar_y
            a_p[y, x - 1] = d[y, x - 1] * detar_y
            a_p[y + 1, x] = d[y + 1, x] * detar_x
            a_p[y - 1, x] = d[y - 1, x] * detar_x

            a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
            b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

            p_add[y, x] = (
                              a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1]
                              + a_p[y + 1, x] * p_add[y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add
                          ) / a_p[y, x]

            p[y, x] = pp[y, x] + alpha_p * p_add[y, x]

    "------------------------------------------u速度修正---------------------------------------------"
    for y in range(n):
        for x in range(n):
            '左下角点'
            if x == 1 and y == 0:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

            '左上角点'
            if x == 1 and y == n - 1:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

            '右下角点'
            if x == n - 1 and y == 0:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

            '右上角点'
            if x == n - 1 and y == n - 1:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * 1.5 * detar_x / detar_y

            '上边界点'
            if y == n - 1 and x != 0 and x != 1 and x != n - 1 and x != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * detar_x / detar_y

            '下边界点'
            if y == 0 and x != 0 and x != 1 and x != n - 1 and x != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 3 / Re * detar_x / detar_y

            '右边界点'
            if x == n - 1 and y != n - 1 and y != 0 and y != n:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max(uu[y, x + 1], 0) + max(-(uu[y, x - 1] + uu[y, x]) / 2, 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 1.5 * detar_x / detar_y

            '左边界点'
            if x == 1 and y != n - 1 and y != 0 and y != n:
                a_p[y, x] = 1.5 * detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-uu[y, x - 1], 0)) \
                            + 1.5 * detar_x * (
                max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 1.5 * detar_x / detar_y

            '内点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != n - 1 and y != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (max((uu[y, x] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y, x - 1]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y + 1, x - 1]) / 2, 0) + max(-(vv[y, x] + vv[y, x - 1]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            d[y, x - 1] = detar_y * alpha_u / a_p[y, x]

            u[y, x] = u[y, x] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])

    "------------------------------------------v速度修正---------------------------------------------"
    for y in range(n):
        for x in range(n):
            '左下角点'
            if x == 0 and y == 1:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-vv[y - 1, x], 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            '左上角点'
            if x == 0 and y == n - 1:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            '右下角点'
            if x == n - 1 and y == 1:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            '右上角点'
            if x == n - 1 and y == n - 1:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            '左边界'
            if x == 0 and y != 0 and y != 1 and y != n - 1 and y != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            '右边界'
            if x == n - 1 and y != 0 and y != 1 and y != n - 1 and y != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 3 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            '下边界'
            if y == 1 and x != 0 and x != n - 1 and x != n:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-vv[y - 1, x], 0)) \
                            + 2 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            '上边界'
            if y == n - 1 and x != 0 and x != n - 1 and x != n:
                a_p[y, x] = detar_x * 1.5 * detar_y / detar_t \
                            + 1.5 * detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max(vv[y + 1, x], 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 2 / Re * 1.5 * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            '内节点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != 1 and y != n - 1 and y != n:
                a_p[y, x] = detar_x * detar_y / detar_t \
                            + detar_y * (
                max((uu[y - 1, x + 1] + uu[y, x + 1]) / 2, 0) + max(-(uu[y, x] + uu[y - 1, x]) / 2, 0)) \
                            + detar_x * (max((vv[y + 1, x] + vv[y, x]) / 2, 0) + max(-(vv[y, x] + vv[y - 1, x]) / 2, 0)) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * detar_x / detar_y

            d[y - 1, x] = detar_x * alpha_v / a_p[y, x]

            v[y, x] = v[y, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

    uu = u
    vv = v
    pp = p
    print(u)

"------------------------------------------写入数据---------------------------------------------"
print(m)
print(u)

for y in range(n):
    for x in range(n+1):
        show_u[m, 0] = x*detar_x
        show_u[m, 1] = y*detar_y
        show_u[m, 2] = u[y, x]
        m += 1
        # print(k)

m = 0
data = open('u.dat','w+')
print('TITLE="u FUNCTION"', file=data)
print('VARIABLES="X","Y","u"', file=data)
print('ZONE T="u" i=%f j=%f c=black' % (n+1, n), file=data)
for m in range((n+1)*(n+1)):
    show_u[m, 0] = round(show_u[m, 0], 3)
    float(show_u[m, 0])
    show_u[m, 1] = round(show_u[m, 1], 3)
    float(show_u[m, 1])
    show_u[m, 2] = round(show_u[m, 2], 3)
    float(show_u[m, 2])
    print('%f     %f     %f' % (show_u[m,0], show_u[m, 1], show_u[m, 2]), file=data)

data.close()