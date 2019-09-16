import numpy as np
"---------------------------------------------Problem-------------------------------------------------------------------"

'假设0＜=x,y<=1的方腔内充满粘性不可压缩流体，左右，下壁固定，上壁以u=1运动，试求Re=100,200,400时的定常解'

"------------------------------------------Initialization--------------------------------------------------------------"
l = 1  # x,y的长度
n = 10  # 内节点法,x,y方向各有１０个节点
detar_x = l / n
detar_y = l / n
detar_t = 1
Re = 100

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
p_add = np.zeros((n+1, n+1))
a_p = np.zeros((n+1, n+1))

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
alpha_p = 0.3
"------------------------------------------Iterative_Solution----------------------------------------------------------"
u[n] = 1


def u_solution (x, y, i, j, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u):
    '设置detar_x与detar_y的长度，主要用于边界点和角点，当为左右边界１．５detar_x，当为上下边界1.5detar_y'
    detar_x = i * detar_x
    detar_y = j * detar_y

    # u_xmax = (max(uu[y, x], 0) - max(-uu[y, x + 1], 0)) - (max(uu[y, x - 1], 0) - max(-uu[y, x], 0))
    # v_ymax = (max(vv[y, x], 0) - max(-vv[y + 1, x], 0)) - (max(vv[y - 1, x], 0) - max(-vv[y, x], 0))

    a_p[y, x] = detar_x * detar_y / detar_t \
                + detar_y * (max(uu[y, x], 0) + max(-uu[y, x], 0)) \
                + detar_x * (max(vv[y, x], 0) + max(-vv[y, x], 0)) \
                + 2 / Re * detar_y / detar_x \
                + 2 / Re * detar_x / detar_y

    a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-uu[y, x+1], 0)
    a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-vv[y+1, x], 0)
    a_w[y, x] = -detar_y / detar_x / Re - detar_y * max(uu[x-1, y], 0)
    a_s[y, x] = -detar_x / detar_y / Re - detar_x * max(vv[y-1, x], 0)

    b = -detar_x * detar_y / detar_t * uu[y, x]
    P = detar_y * (p[y, x] - p[y, x - 1])

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    a_e[y, x] = a_e[y, x] * E
    a_n[y, x] = a_n[y, x] * N
    a_w[y, x] = a_w[y, x] * W
    a_s[y, x] = a_s[y, x] * S

    if N == 0:
        b += -detar_x / detar_y / Re - detar_x * max(-vv[y+1, x], 0)

    '亚松弛设置,注意两者的顺序'
    b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
    a_p[y, x] = a_p[y, x] / alpha_u

    u[y, x] = (a_e[y, x] * u[y, x + 1]
               + a_n[y, x] * u[y + 1, x]
               + a_w[y, x] * u[y, x - 1]
               + a_s[y, x] * u[y - 1, x]
               + b + P) / (-a_p[y, x])

    return u[y,x]

def v_solution (x, y, i, j, E, W, N, S, detar_x, detar_y, Re, vv, uu, alpha_v):
    '设置detar_x与detar_y的长度，主要用于边界点和角点，当为左右边界１．５detar_x，当为上下边界1.5detar_y'
    detar_x = i * detar_x
    detar_y = j * detar_y

    # u_xmax = (max(uu[y, x], 0) - max(-uu[y, x + 1], 0)) - (max(uu[y, x - 1], 0) - max(-uu[y, x], 0))
    # v_ymax = (max(vv[y, x], 0) - max(-vv[y + 1, x], 0)) - (max(vv[y - 1, x], 0) - max(-vv[y, x], 0))

    a_p[y, x] = detar_x * detar_y / detar_t \
                + detar_x * (max(vv[y, x], 0) + max(-vv[y, x], 0)) \
                + detar_y * (max(uu[y, x], 0) + max(-uu[y, x], 0)) \
                + 2 / Re * detar_y / detar_x \
                + 2 / Re * detar_x / detar_y

    a_e[y, x] = -detar_y / detar_x / Re - detar_y * max(-uu[y, x+1], 0)
    a_n[y, x] = -detar_x / detar_y / Re - detar_x * max(-vv[y+1, x], 0)
    a_w[y, x] = -detar_y / detar_x / Re - detar_y * max(uu[y, x-1], 0)
    a_s[y, x] = -detar_x / detar_y / Re - detar_x * max(vv[y-1, x], 0)

    b = -detar_x * detar_y / detar_t * vv[y, x]
    P = detar_x * (p[y, x] - p[y - 1, x])

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    a_e[y, x] = a_e[y, x] * E
    a_n[y, x] = a_n[y, x] * N
    a_w[y, x] = a_w[y, x] * W
    a_s[y, x] = a_s[y, x] * S

    '亚松弛设置'
    b = b + (1 - alpha_v) * a_p[y, x] / alpha_v * vv[y, x]
    a_p[y, x] = a_p[y, x] / alpha_v

    v[y, x] = (a_e[y, x] * v[y, x + 1]
               + a_n[y, x] * v[y + 1, x]
               + a_w[y, x] * v[y, x - 1]
               + a_s[y, x] * v[y - 1, x]
               + b + P) / (-a_p[y, x])

    return v[y,x]

def p_add_solution (x, y, i_u, j_u, i_v, j_v, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p):
    detar_x_u = i_u * detar_x
    detar_y_u = j_u * detar_y

    detar_x_v = i_v * detar_x
    detar_y_v = j_v * detar_y

    # u_xmax_W = (max(uu[y, x], 0) - max(-uu[y, x + 1], 0)) - (max(uu[y, x - 1], 0) - max(-uu[y, x], 0))
    # v_ymax_W = (max(vv[y, x], 0) - max(-vv[y + 1, x], 0)) - (max(vv[y - 1, x], 0) - max(-vv[y, x], 0))
    #
    # u_xmax_S = (max(uu[y, x], 0) - max(-uu[y, x + 1], 0)) - (max(uu[y, x - 1], 0) - max(-uu[y, x], 0))
    # v_ymax_S = (max(vv[y, x], 0) - max(-vv[y + 1, x], 0)) - (max(vv[y - 1, x], 0) - max(-vv[y, x], 0))
    #
    # if x == n-1:
    #     u_xmax_E = - (max(uu[y, x], 0) - max(-uu[y, x+1], 0))
    # if x != n-1:
    #     u_xmax_E = (max(uu[y, x+1], 0) - max(-uu[y, x + 2], 0)) - (max(uu[y, x], 0) - max(-uu[y, x+1], 0))
    # v_ymax_E = (max(vv[y, x+1], 0) - max(-vv[y + 1, x+1], 0)) - (max(vv[y - 1, x+1], 0) - max(-vv[y, x+1], 0))
    #
    # u_xmax_N = (max(uu[y+1, x], 0) - max(-uu[y+1, x + 1], 0)) - (max(uu[y+1, x-1], 0) - max(-uu[y+1, x], 0))
    # if y == n-1:
    #     v_ymax_N = - (max(vv[y, x], 0) - max(-vv[y+1, x], 0))
    # if y != n-1:
    #     v_ymax_N = (max(vv[y + 1, x], 0) - max(-vv[y + 2, x], 0)) - (max(vv[y, x], 0) - max(-vv[y + 1, x], 0))

    if x == n - 1:
        a_e[y, x] = detar_x_u * detar_y_u / detar_t \
                    + 2 / Re * detar_y_u / detar_x_u \
                    + 2 / Re * detar_x_u / detar_y_u
    else:
        a_e[y, x] = detar_x_u * detar_y_u / detar_t \
                    + detar_y_u * (max(uu[y, x+1], 0) + max(-uu[y, x+1], 0)) \
                    + detar_x_u * (max(vv[y, x+1], 0) + max(-vv[y, x+1], 0)) \
                    + 2 / Re * detar_y_u / detar_x_u \
                    + 2 / Re * detar_x_u / detar_y_u

    a_w[y, x] = detar_x_u * detar_y_u / detar_t \
                + detar_y_u * (max(uu[y, x], 0) + max(-uu[y, x], 0)) \
                + detar_x_u * (max(vv[y, x], 0) + max(-vv[y, x], 0)) \
                + 2 / Re * detar_y_u / detar_x_u \
                + 2 / Re * detar_x_u / detar_y_u

    if y == n-1:
        a_n[y, x] = detar_x_v * detar_y_v / detar_t \
                    + 2 / Re * detar_y_v / detar_x_v \
                    + 2 / Re * detar_x_v / detar_y_v
    else:
        a_n[y, x] = detar_x_v * detar_y_v / detar_t \
                    + detar_x_v * (max(vv[y+1, x], 0) + max(-vv[y+1, x], 0)) \
                    + detar_y_v * (max(uu[y+1, x], 0) + max(-uu[y+1, x], 0)) \
                    + 2 / Re * detar_y_v / detar_x_v \
                    + 2 / Re * detar_x_v / detar_y_v

    a_s[y, x] = detar_x_v * detar_y_v / detar_t \
                + detar_x_v * (max(vv[y, x], 0) + max(-vv[y, x], 0)) \
                + detar_y_v * (max(uu[y, x], 0) + max(-uu[y, x], 0)) \
                + 2 / Re * detar_y_v / detar_x_v \
                + 2 / Re * detar_x_v / detar_y_v

    d[y, x + 1] = detar_y_u * alpha_u / a_e[y, x]
    d[y, x - 1] = detar_y_u * alpha_u / a_w[y, x]
    d[y + 1, x] = detar_x_v * alpha_v / a_n[y, x]
    d[y - 1, x] = detar_x_v * alpha_v / a_s[y, x]

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    d[y, x + 1] = d[y, x + 1] * E
    d[y + 1, x] = d[y + 1, x] * N
    d[y, x - 1] = d[y, x - 1] * W
    d[y - 1, x] = d[y - 1, x] * S

    a_p[y, x + 1] = d[y, x + 1] * detar_y_u
    a_p[y, x - 1] = d[y, x - 1] * detar_y_u
    a_p[y+1, x] = d[y+1, x] * detar_x_v
    a_p[y-1, x] = d[y-1, x] * detar_x_v

    a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
    b_add = (u[y, x] - u[y, x + 1]) * detar_y_u + (v[y, x] - v[y + 1, x]) * detar_x_v

    p_add[y, x] = (
                    a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1]
                    + a_p[y + 1, x] *p_add[y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add
                  ) / a_p[y, x]

    p[y, x] = p[y, x] + alpha_p * p_add[y, x]

    return p[y,x]

def u_add_solution (x, y, i_u, j_u, i_v, j_v, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p):
    detar_x_u = i_u * detar_x
    detar_y_u = j_u * detar_y

    # detar_x_v = i_v * detar_x
    # detar_y_v = j_v * detar_y

    # u_xmax_W = (max(uu[y, x], 0) - max(-uu[y, x + 1], 0)) - (max(uu[y, x - 1], 0) - max(-uu[y, x], 0))
    # v_ymax_W = (max(vv[y, x], 0) - max(-vv[y + 1, x], 0)) - (max(vv[y - 1, x], 0) - max(-vv[y, x], 0))

    # u_xmax_S = (max(uu[y, x], 0) - max(-uu[y, x + 1], 0)) - (max(uu[y, x - 1], 0) - max(-uu[y, x], 0))
    # v_ymax_S = (max(vv[y, x], 0) - max(-vv[y + 1, x], 0)) - (max(vv[y - 1, x], 0) - max(-vv[y, x], 0))

    # if x == n-1:
    #     u_xmax_E = - (max(uu[y, x], 0) - max(-uu[y, x+1], 0))
    # if x != n-1:
    #     u_xmax_E = (max(uu[y, x+1], 0) - max(-uu[y, x + 2], 0)) - (max(uu[y, x], 0) - max(-uu[y, x+1], 0))
    # v_ymax_E = (max(vv[y, x+1], 0) - max(-vv[y + 1, x+1], 0)) - (max(vv[y - 1, x+1], 0) - max(-vv[y, x+1], 0))
    #
    # u_xmax_N = (max(uu[y+1, x], 0) - max(-uu[y+1, x + 1], 0)) - (max(uu[y+1, x-1], 0) - max(-uu[y+1, x], 0))
    # if y == n-1:
    #     v_ymax_N = - (max(vv[y, x], 0) - max(-vv[y+1, x], 0))
    # if y != n-1:
    #     v_ymax_N = (max(vv[y + 1, x], 0) - max(-vv[y + 2, x], 0)) - (max(vv[y, x], 0) - max(-vv[y + 1, x], 0))

    # a_e[y, x] = detar_x_u * detar_y_u / detar_t \
    #                     + detar_y_u * u_xmax_E \
    #                     + 2 / Re * detar_y_u / detar_x_u \
    #                     + 2 / Re * detar_x_u / detar_y_u \
    #                     + v_ymax_E * detar_x_u

    a_w[y, x] = detar_x_u * detar_y_u / detar_t \
                + detar_y_u * (max(uu[y, x], 0) + max(-uu[y, x], 0)) \
                + detar_x_u * (max(vv[y, x], 0) + max(-vv[y, x], 0)) \
                + 2 / Re * detar_y_u / detar_x_u \
                + 2 / Re * detar_x_u / detar_y_u

    # a_n[y, x] = detar_x_v * detar_y_v / detar_t \
    #                 + detar_x_v * v_ymax_N \
    #                 + 2 / Re * detar_y_v / detar_x_v \
    #                 + 2 / Re * detar_x_v / detar_y_v \
    #                 + detar_y_v * u_xmax_N

    # a_s[y, x] = detar_x_v * detar_y_v / detar_t \
    #                 + detar_x_v * v_ymax_S \
    #                 + 2 / Re * detar_y_v / detar_x_v \
    #                 + 2 / Re * detar_x_v / detar_y_v \
    #                 + detar_y_v * u_xmax_S


    # d[y, x + 1] = detar_y_u * alpha_p / a_e[y, x]
    d[y, x - 1] = detar_y_u * alpha_u / a_w[y, x]
    # d[y + 1, x] = detar_x_v * alpha_p / a_n[y, x]
    # d[y - 1, x] = detar_x_v * alpha_p / a_s[y, x]

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    # d[y, x + 1] = d[y, x + 1] * E
    # d[y + 1, x] = d[y + 1, x] * N
    d[y, x - 1] = d[y, x - 1] * W
    # d[y - 1, x] = d[y - 1, x] * S

    u[y, x] = u[y, x] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
    # v[y, x] = v[y, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

    return u[y, x]

def v_add_solution (x, y, i_u, j_u, i_v, j_v, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p):
    # detar_x_u = i_u * detar_x
    # detar_y_u = j_u * detar_y

    detar_x_v = i_v * detar_x
    detar_y_v = j_v * detar_y

    # u_xmax_W = (max(uu[y, x], 0) - max(-uu[y, x + 1], 0)) - (max(uu[y, x - 1], 0) - max(-uu[y, x], 0))
    # v_ymax_W = (max(vv[y, x], 0) - max(-vv[y + 1, x], 0)) - (max(vv[y - 1, x], 0) - max(-vv[y, x], 0))

    u_xmax_S = (max(uu[y, x], 0) - max(-uu[y, x + 1], 0)) - (max(uu[y, x - 1], 0) - max(-uu[y, x], 0))
    v_ymax_S = (max(vv[y, x], 0) - max(-vv[y + 1, x], 0)) - (max(vv[y - 1, x], 0) - max(-vv[y, x], 0))

    # if x == n-1:
    #     u_xmax_E = - (max(uu[y, x], 0) - max(-uu[y, x+1], 0))
    # if x != n-1:
    #     u_xmax_E = (max(uu[y, x+1], 0) - max(-uu[y, x + 2], 0)) - (max(uu[y, x], 0) - max(-uu[y, x+1], 0))
    # v_ymax_E = (max(vv[y, x+1], 0) - max(-vv[y + 1, x+1], 0)) - (max(vv[y - 1, x+1], 0) - max(-vv[y, x+1], 0))
    #
    # u_xmax_N = (max(uu[y+1, x], 0) - max(-uu[y+1, x + 1], 0)) - (max(uu[y+1, x-1], 0) - max(-uu[y+1, x], 0))
    # if y == n-1:
    #     v_ymax_N = - (max(vv[y, x], 0) - max(-vv[y+1, x], 0))
    # if y != n-1:
    #     v_ymax_N = (max(vv[y + 1, x], 0) - max(-vv[y + 2, x], 0)) - (max(vv[y, x], 0) - max(-vv[y + 1, x], 0))

    # a_e[y, x] = detar_x_u * detar_y_u / detar_t \
    #                     + detar_y_u * u_xmax_E \
    #                     + 2 / Re * detar_y_u / detar_x_u \
    #                     + 2 / Re * detar_x_u / detar_y_u \
    #                     + v_ymax_E * detar_x_u

    # a_w[y, x] = detar_x_u * detar_y_u / detar_t \
    #                 + detar_y_u * u_xmax_W \
    #                 + 2 / Re * detar_y_u / detar_x_u \
    #                 + 2 / Re * detar_x_u / detar_y_u \
    #                 + v_ymax_W * detar_x_u

    # a_n[y, x] = detar_x_v * detar_y_v / detar_t \
    #                 + detar_x_v * v_ymax_N \
    #                 + 2 / Re * detar_y_v / detar_x_v \
    #                 + 2 / Re * detar_x_v / detar_y_v \
    #                 + detar_y_v * u_xmax_N

    a_s[y, x] = detar_x_v * detar_y_v / detar_t \
                + detar_x_v * (max(vv[y, x], 0) + max(-vv[y, x], 0)) \
                + detar_y_v * (max(uu[y, x], 0) + max(-uu[y, x], 0)) \
                + 2 / Re * detar_y_v / detar_x_v \
                + 2 / Re * detar_x_v / detar_y_v


    # d[y, x + 1] = detar_y_u * alpha_p / a_e[y, x]
    # d[y, x - 1] = detar_y_u * alpha_p / a_w[y, x]
    # d[y + 1, x] = detar_x_v * alpha_p / a_n[y, x]
    d[y - 1, x] = detar_x_v * alpha_v / a_s[y, x]

    '设置ue,un,uw,ns的系数，主要用于边界点和角点，左边界a_w=0，右边界a_e=0,上边界a_n=0,下边界a_s=0'
    # d[y, x + 1] = d[y, x + 1] * E
    # d[y + 1, x] = d[y + 1, x] * N
    # d[y, x - 1] = d[y, x - 1] * W
    d[y - 1, x] = d[y - 1, x] * S

    # u[y, x] = u[y, x] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
    v[y, x] = v[y, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

    return v[y, x]

for k in range(1000):
    "------------------------------------------u控制方程计算---------------------------------------------"
    for y in range(n+1):
        for x in range(n+1):

            '四个脚点'
            if x == 1 and y == 0:
                u_left_down = u_solution(x, y, 1.5, j, E, 0, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u)
                u[y, x] = u_left_down

            if x == 1 and y == n - 1:
                u_left_up = u_solution(x, y, 1.5, j, E, 0, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u)
                u[y, x] = u_left_up

            if x == n - 1 and y == 0:
                u_right_down = u_solution(x, y, 1.5, j, 0, W, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u)
                u[y, x] = u_right_down

            if x == n - 1 and y == n - 1:
                u_right_up = u_solution(x, y, 1.5, j, 0, W, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u)
                u[y, x] = u_right_up

            '四个边界'
            if x == 1 and y != n - 1 and y != 0 and y != n:
                u_left = u_solution(x, y, 1.5, j, E, 0, N, S, detar_x, detar_y, Re, uu, vv, alpha_u)
                u[y, x] = u_left

            if x == n - 1 and y != n - 1 and y != 0 and y != n:
                u_right = u_solution(x, y, 1.5, j, 0, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u)
                u[y, x] = u_right

            if y == 0 and x != 0 and x != 1 and x != n - 1 and x != n:
                u_down = u_solution(x, y, i, j, E, W, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u)
                u[y, x] = u_down

            if y == n - 1 and x != 0 and x != 1 and x != n - 1 and x != n:
                u_up = u_solution(x, y, i, j, E, W, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u)
                u[y, x] = u_up

            '内节点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != n - 1 and y != n:
                u_in = u_solution(x, y, i, j, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u)
                u[y, x] = u_in

    "------------------------------------------v控制方程计算---------------------------------------------"
    for y in range(n+1):
        for x in range(n+1):

            '四个角点'
            if x == 0 and y == 1:
                v_left_down = v_solution(x, y, i, 1.5, E, 0, N, 0, detar_x, detar_y, Re, vv, uu, alpha_v)
                v[y, x] = v_left_down

            if x == 0 and y == n - 1:
                v_left_up = v_solution(x, y, i, 1.5, E, 0, 0, S, detar_x, detar_y, Re, vv, uu, alpha_v)
                v[y, x] = v_left_up

            if x == n - 1 and y == 1:
                v_right_down = v_solution(x, y, i, 1.5, 0, W, N, 0, detar_x, detar_y, Re, vv, uu, alpha_v)
                v[y, x] = v_right_down

            if x == n - 1 and y == n - 1:
                v_right_up = v_solution(x, y, i, 1.5, 0, W, 0, S, detar_x, detar_y, Re, vv, uu, alpha_v)
                v[y, x] = v_right_up

            '四个边界'
            if x == 0 and y != 0 and y != 1 and y != n - 1 and y != n:
                v_left = v_solution(x, y, i, j, E, 0, N, S, detar_x, detar_y, Re, vv, uu, alpha_v)
                v[y, x] = v_left

            if x == n - 1 and y != 0 and y != 1 and y != n - 1 and y != n:
                v_right = v_solution(x, y, i, j, 0, W, N, S, detar_x, detar_y, Re, vv, uu, alpha_v)
                v[y, x] = v_right

            if y == 1 and x != 0 and x != n - 1 and x != n:
                v_down = v_solution(x, y, i, 1.5, E, W, N, 0, detar_x, detar_y, Re, vv, uu, alpha_v)
                v[y, x] = v_down

            if y == n - 1 and x != 0 and x != n - 1 and x != n:
                v_up = v_solution(x, y, i, 1.5, E, W, 0, S, detar_x, detar_y, Re, vv, uu, alpha_v)
                v[y, x] = v_up

            '内节点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != 1 and y != n - 1 and y != n:
                v_in = v_solution(x, y, i, j, E, W, N, S, detar_x, detar_y, Re, vv, uu, alpha_v)
                v[y, x] = v_in

    "------------------------------------------压力修正---------------------------------------------"
    for y in range(n+1):
        for x in range(n+1):
            '四个角点'
            if x == 0 and y == 0:
                p_left_down = p_add_solution(x, y, 1.5, j_u, i_v, 1.5, E, 0, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                p[y, x] = p_left_down

            if x == n - 1 and y == 0:
                p_right_down = p_add_solution(x, y, i_u, j_u, i_v, j_v, 0, W, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                p[y, x] = p_right_down

            if x == 0 and y == n - 1:
                p_left_up = p_add_solution(x, y, 1.5, j_u, i_v, 1.5, E, 0, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                p[y, x] = p_left_up

            if x == n - 1 and y == n - 1:
                p_right_up = p_add_solution(x, y, 1.5, j_u, i_v, 1.5, 0, W, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                p[y, x] = p_right_up

            '四个边界'
            if x == 0 and y != 0 and y != n - 1 and y != n:
                p_left = p_add_solution(x, y, 1.5, j_u, i_v, j_v, E, 0, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                p[y, x] = p_left
            if x == n - 1 and y != 0 and y != n - 1 and y != n:
                p_right = p_add_solution(x, y, 1.5, j_u, i_v, j_v, 0, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                p[y, x] = p_right
            if y == 0 and x != 0 and x != n - 1 and x != n:
                p_down = p_add_solution(x, y, i_u, j_u, i_v, 1.5, E, W, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                p[y, x] = p_down
            if y == n - 1 and x != 0 and x != n - 1 and x != n:
                p_up = p_add_solution(x, y, i_u, j_u, i_v, 1.5, E, W, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                p[y, x] = p_up

            '内部节点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != n - 1 and y != n:
                p_in = p_add_solution(x, y, i_u, j_u, i_v, j_v, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                p[y, x] = p_in

    "------------------------------------------速度修正---------------------------------------------"
    for y in range(n+1):
        for x in range(n+1):
            '四个角点'
            if x == 0 and y == 0:
                # u_add_left_down = u_add_solution(x, y, 1.5, j_u, i_v, 1.5, E, 0, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                u[y, x] = 0

                # v_add_left_down = v_add_solution(x, y, 1.5, j_u, i_v, 1.5, E, 0, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                v[y, x] = 0

            if x == n - 1 and y == 0:
                u_add_right_down = u_add_solution(x, y, 1.5, j_u, i_v, 1.5, 0, W, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                u[y, x] = u_add_right_down

                # v_add_right_down = v_add_solution(x, y, 1.5, j_u, i_v, 1.5, 0, W, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                v[y, x] = 0

            if x == 0 and y == n - 1:
                # u_add_left_top = u_add_solution(x, y, 1.5, j_u, i_v, 1.5, E, 0, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                u[y, x] = 0

                v_add_left_top = v_add_solution(x, y, 1.5, j_u, i_v, 1.5, E, 0, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                v[y, x] = v_add_left_top

            if x == n - 1 and y == n - 1:
                u_add_right_top = u_add_solution(x, y, 1.5, j_u, i_v, 1.5, 0, W, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                u[y, x] = u_add_right_top

                v_add_right_top = v_add_solution(x, y, 1.5, j_u, i_v, 1.5, 0, W, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                v[y, x] = v_add_right_top

            '四个边界'
            if x == 0 and y != 0 and y != n - 1 and y != n:
                # u_add_left = u_add_solution(x, y, 1.5, j_u, i_v, j_v, E, 0, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                u[y, x] = 0

                v_add_left = v_add_solution(x, y, 1.5, j_u, i_v, j_v, E, 0, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                v[y, x] = v_add_left

            if x == n - 1 and y != 0 and y != n - 1 and y != n:
                u_add_right = u_add_solution(x, y, 1.5, j_u, i_v, j_v, 0, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                u[y, x] = u_add_right

                v_add_right = v_add_solution(x, y, 1.5, j_u, i_v, j_v, 0, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                v[y, x] = v_add_right

            if y == 0 and x != 0 and x != n - 1 and x != n:
                u_add_down = u_add_solution(x, y, i_u, j_u, i_v, 1.5, E, W, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                u[y, x] = u_add_down

                # v_add_down = v_add_solution(x, y, i_u, j_u, i_v, 1.5, E, W, N, 0, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                v[y, x] = 0

            if y == n - 1 and x != 0 and x != n - 1 and x != n:
                u_add_top = u_add_solution(x, y, i_u, j_u, i_v, 1.5, E, W, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                u[y, x] = u_add_top

                v_add_top = v_add_solution(x, y, i_u, j_u, i_v, 1.5, E, W, 0, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                v[y, x] = v_add_top

            '内节点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != n - 1 and y != n:
                u_add_in = u_add_solution (x, y, i_u, j_u, i_v, j_v, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                u[y, x] = u_add_in

                v_add_in = v_add_solution (x, y, i_u, j_u, i_v, j_v, E, W, N, S, detar_x, detar_y, Re, uu, vv, alpha_u, alpha_v, alpha_p)
                # v[y, x] = v_add_in

    uu = u
    vv = v
    print(u)
    # print(v)
    # print(p)