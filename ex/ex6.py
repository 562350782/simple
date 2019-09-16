import numpy as np

"---------------------------------------------Problem-------------------------------------------------------------------"

'假设0＜=x,y<=1的方腔内充满粘性不可压缩流体，左右，下壁固定，上壁以u=1运动，试求Re=100,200,400时的定常解'

"------------------------------------------Initialization--------------------------------------------------------------"
l = 1  # x,y的长度
n = 5  # 内节点法,x,y方向各有１０个节点
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

a_v = np.zeros((n+1, n+1))  # n+1时层v的各项系数
aa_v = np.zeros((n+1, n+1))  # n时层v的各项系数
a_v_add = np.zeros((n+1, n+1))

v = np.zeros((n+1, n+1))  # n时层次的速度u
vv = np.zeros((n+1, n+1))  # n+1时层次的速度u
v_add = np.zeros((n+1, n+1))

p = np.zeros((n+1, n+1))
p_add = np.zeros((n+1, n+1))
a_p = np.zeros((n+1, n+1))

d = np.zeros((n+1, n+1))

'亚松弛'
alpha_u = 0.5
alpha_v = 0.5
alpha_p = 0.3
"------------------------------------------Iterative_Solution----------------------------------------------------------"

u[5] = 1
print(u)

for k in range(100):
    for y in range(n+1):
        for x in range(n+1):
            "------------------------------------------u控制方程计算---------------------------------------------"
            '四个脚点'
            if x == 1 and y == 0:
                a_u[y, x] = 3 / 2 * detar_x * detar_y / detar_t \
                            + uu[y, x] * detar_y \
                            + 3 * detar_x / Re / detar_y \
                            + 3 / 8 * detar_x * (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) \
                            + 2 / Re * detar_y / detar_x

                a_u[y, x + 1] = - detar_y / detar_x / Re
                a_u[y, x - 1] = 0
                a_u[y + 1, x] = - 3 / 2 * detar_x / detar_y / Re
                a_u[y - 1, x] = 0

                b = - 3 / 2 * uu[y, x] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
                a_u[y, x] = -a_u[y, x] / alpha_u

                u[y, x] = (
                              a_u[y, x + 1] * u[y, x + 1]
                              + a_u[y, x - 1] * u[y, x - 1]
                              + a_u[y + 1, x] * u[y + 1, x]
                              + a_u[y - 1, x] * u[y - 1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y, x-1]) * detar_y
                          ) / a_u[y, x]

            if x == 1 and y == n-1:
                a_u[y, x] = 3 / 2 * detar_x * detar_y / detar_t \
                            + uu[y, x] * detar_y \
                            + 3 * detar_x / Re / detar_y \
                            + 3 / 8 * detar_x * (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) \
                            + 2 / Re * detar_y / detar_x

                a_u[x, y+1] = - detar_y / detar_x / Re
                a_u[x, y-1] = 0
                a_u[y+1, x] = 0
                a_u[y-1, x] = - 3 / 8 * detar_x * (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) \
                                - 3 / 2 / Re * detar_x / detar_y

                b = - 3 / 2 * uu[y, x] * detar_x * detar_y / detar_t - 3/2*detar_x/Re/detar_y
                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
                a_u[y, x] = -a_u[y, x] / alpha_u

                u[y, x] = (
                              a_u[x, y+1] * u[x, y+1]
                              + a_u[x, y-1] * u[x, y-1]
                              + a_u[y+1, x] * u[y+1, x]
                              + a_u[y-1, x] * u[y-1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y, x-1]) * detar_y
                          ) / a_u[y, x]

            if x == n-1 and y == 0:
                a_u[y, x] = 3 / 2 * detar_x * detar_y / detar_t \
                            + 3 / 8 * detar_x * (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 3 / 2 * detar_x / detar_y

                a_u[y, x+1] = 0
                a_u[y, x-1] = -uu[y, x] * detar_y - detar_y / Re / detar_x
                a_u[y+1, x] = -3 / 2 * detar_x / detar_y / Re
                a_u[y-1, x] = 0

                b = - 3 / 2 * uu[y, x] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
                a_u[y, x] = -a_u[y, x] / alpha_u

                u[y, x] = (
                              a_u[y, x+1] * u[y, x+1]
                              + a_u[y, x-1] * u[y, x-1]
                              + a_u[y+1, x] * u[y+1, x]
                              + a_u[y-1, x] * u[y-1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y, x-1]) * detar_y
                          ) / a_u[y, x]

            if x == n-1 and y == n-1:
                a_u[y, x] = 3 / 2 * detar_x * detar_y / detar_t \
                            + 3 / 8 * detar_x * (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 3 / 2 * detar_x / detar_y

                a_u[y, x+1] = 0
                a_u[y, x - 1] = -uu[y, x] * detar_y - detar_y / Re / detar_x
                a_u[y + 1, x] = 0
                a_u[y - 1, x] = -3 / 8 * detar_x * (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) \
                                - 3 / 2 * detar_x / Re / detar_y

                b = - 3 / 2 * uu[y, x] * detar_x * detar_y / detar_t - 3 / 2 * detar_x / detar_y / Re

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
                a_u[y, x] = -a_u[y, x] / alpha_u

                u[y, x] = (
                              a_u[y, x+1] * u[y, x+1]
                              + a_u[y, x - 1] * u[y, x - 1]
                              + a_u[y + 1, x] * u[y + 1, x]
                              + a_u[y - 1, x] * u[y - 1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y, x-1]) * detar_y
                        ) / a_u[y, x]

            '四个边界'
            if x == 1 and y != n - 1 and y != 0 and y != n:
                a_u[y, x] = 3 / 2 * detar_x * detar_y / detar_t \
                            + uu[y, x] * detar_y \
                            + 3 * detar_x / Re / detar_y \
                            + 3 / 8 * detar_x * (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) \
                            + 2 / Re * detar_y / detar_x

                a_u[y, x+1] = - detar_y / detar_x / Re
                a_u[y, x - 1] = 0
                a_u[y + 1, x] = - 3 / 2 * detar_x / detar_y / Re
                a_u[y - 1, x] = - 3 / 8 * detar_x * (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) \
                                - 3 / 2 / Re * detar_x / detar_y

                b = - 3 / 2 * uu[y, x] * detar_x * detar_y / detar_t

                # '亚松弛设置,注意两者的顺序'
                # b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
                # a_u[y, x] = -a_u[y, x] / alpha_u

                u[y, x] = (
                              a_u[y, x+1] * u[y, x+1]
                              + a_u[y, x - 1] * u[y, x - 1]
                              + a_u[y + 1, x] * u[y + 1, x]
                              + a_u[y - 1, x] * u[y - 1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y, x - 1]) * detar_y
                          ) / a_u[y, x]

            if x == n-1 and y!= n-1 and y != 0 and y!=n:
                a_u[y, x] = 3 / 2 * detar_x * detar_y / detar_t \
                            + 3 / 8 * detar_x * (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 3 / 2 * detar_x / detar_y

                a_u[y, x+1] = 0
                a_u[y, x - 1] = -uu[y, x] * detar_y - detar_y / Re / detar_x
                a_u[y+1, x] = -3 / 2 * detar_x / detar_y / Re
                a_u[y-1, x] = -3 / 8 * detar_x * (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) \
                                - 3 / 2 * detar_x / Re / detar_y

                b = - 3 / 2 * uu[y, x] * detar_x * detar_y / detar_t

                # '亚松弛设置,注意两者的顺序'
                # b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
                # a_u[y, x] = -a_u[y, x] / alpha_u

                u[y, x] = (
                              a_u[y, x+1] * u[y, x+1]
                              + a_u[y, x - 1] * u[y, x - 1]
                              + a_u[y+1, x] * u[y+1, x]
                              + a_u[y-1, x] * u[y-1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y, x-1]) * detar_y
                          ) / a_u[y, x]

            if y == 0 and x!=0 and x!= 1 and x != n-1 and x!=n:
                a_u[y, x] = (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) / 4 * detar_x \
                            + uu[y, x] * detar_y \
                            + detar_x * detar_y / detar_t \
                            + detar_y / detar_x * 2 / Re \
                            + detar_x / detar_y * 2 / Re

                a_u[y, x+1] = - detar_y / detar_x / Re
                a_u[y, x-1] = - uu[y, x] * detar_y \
                                - detar_y / detar_x / Re
                a_u[y+1, x] = - detar_x / detar_y / Re
                a_u[y-1, x] = 0

                b = - uu[y, x] * detar_x * detar_y / detar_t

                # '亚松弛设置,注意两者的顺序'
                # b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
                # a_u[y, x] = -a_u[y, x] / alpha_u

                u[y, x] = (
                              a_u[y, x+1] * u[y, x+1]
                              + a_u[y, x-1] * u[y, x-1]
                              + a_u[y+1, x] * u[y+1, x]
                              + a_u[y-1, x] * u[y-1, x]
                              + b
                              + (p[y, x] - p[y,x-1]) * detar_y
                          ) / a_u[y, x]

            if y == n - 1 and x != 0 and x != 1 and x != n - 1 and x != n:
                a_u[y, x] = (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) / 4 * detar_x \
                            + uu[y, x] * detar_y \
                            + detar_x * detar_y / detar_t \
                            + detar_y / detar_x * 2 / Re \
                            + detar_x / detar_y * 2 / Re

                a_u[y, x + 1] = - detar_y / detar_x / Re
                a_u[y, x - 1] = - uu[y, x] * detar_y \
                                - detar_y / detar_x / Re
                a_u[y + 1, x] = 0
                a_u[y - 1, x] = - (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) / 4 * detar_x \
                                - detar_x / Re / detar_y

                b = - uu[y, x] * detar_x * detar_y / detar_t - 3 / 2 * detar_x / detar_y / Re

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
                a_u[y, x] = -a_u[y, x] / alpha_u

                u[y, x] = (
                              a_u[y, x + 1] * u[y, x + 1]
                              + a_u[y, x - 1] * u[y, x - 1]
                              + a_u[y + 1, x] * u[y + 1, x]
                              + a_u[y - 1, x] * u[y - 1, x]
                              + b
                              + (p[y, x] - p[x, y-1]) * detar_y
                          ) / a_u[y, x]

            '内部节点'
            if x != 0 and x!= 1 and x != n - 1 and x!= n and y != 0 and y != n - 1 and y!= n:
                a_u[y, x] = (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) / 4 * detar_x \
                            + uu[y, x] * detar_y \
                            + detar_x * detar_y / detar_t \
                            + detar_y / detar_x * 2 / Re\
                            + detar_x / detar_y * 2 / Re

                a_u[y, x + 1] = - detar_y / detar_x / Re
                a_u[y, x - 1] = - uu[y, x] * detar_y \
                                - detar_y / detar_x / Re
                a_u[y+1, x] = - detar_x / detar_y / Re
                a_u[y-1, x] = - (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) / 4 * detar_x \
                                - detar_x / Re / detar_y

                b = - uu[y, x] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[y, x] / alpha_u * uu[y, x]
                a_u[y, x] = -a_u[y, x] / alpha_u

                u[y, x] = (
                              a_u[y, x + 1] * u[y, x + 1]
                              + a_u[y, x - 1] * u[y, x - 1]
                              + a_u[y+1, x] * u[y+1, x]
                              + a_u[y-1, x] * u[y-1, x]
                              + b
                              + (p[y, x] - p[y, x - 1]) * detar_y
                          ) / a_u[y, x]
            "------------------------------------------v控制方程计算---------------------------------------------"
            '四个脚点'
            if x == 0 and y == 1:
                a_v[y, x] = vv[y, x] * detar_x \
                            + 3 / 8 * detar_y * (uu[y, x+1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x]) \
                            + 3 / 2 * detar_y * detar_x / detar_t \
                            + 3 * detar_y / detar_x / Re \
                            + 2 * detar_x / detar_y / Re

                a_v[y, x + 1] = - 3 / 2 * detar_y / detar_x / Re
                a_v[y, x - 1] = 0
                a_v[y+1, x] = - detar_x / detar_y / Re
                a_v[y-1, x] = 0

                b = - 3 / 2 * vv[y, x] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[y, x] / alpha_v * vv[y, x]
                a_v[y, x] = -a_v[y, x] / alpha_v

                v[y, x] = (
                              a_v[y, x + 1] * v[y, x + 1]
                              + a_v[y, x - 1] * v[y, x - 1]
                              + a_v[y+1, x] * v[y+1, x]
                              + a_v[y-1, x] * v[y-1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y-1, x]) * detar_x
                          ) / a_v[y, x]

            if x == 0 and y == n-1:
                a_v[y, x] = 3 / 2 * detar_y * detar_x \
                            + vv[y, x] * detar_x * detar_t \
                            + 3 / 8 * detar_y * detar_t * (uu[y, x+1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x]) \
                            + 3 * detar_y * detar_t / detar_x / Re \
                            + 2 * detar_x * detar_t / detar_y / Re

                a_v[y, x+1] = - 3 / 2 * detar_y * detar_t / detar_x / Re
                a_v[y, x-1] = 0
                a_v[y+1, x] = 0
                a_v[y-1, x] = - detar_x * detar_t / detar_y / Re \
                                - vv[y, x] * detar_x * detar_t

                b = - 3 / 2 * vv[y, x] * detar_x * detar_y

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[y, x] / alpha_v * vv[y, x]
                a_v[y, x] = -a_v[y, x] / alpha_v

                v[y, x] = (
                              a_v[y, x+1] * v[y, x+1]
                              + a_v[y, x-1] * v[y, x-1]
                              + a_v[y+1, x] * v[y+1, x]
                              + a_v[y-1, x] * v[y-1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y-1, x]) * detar_x * detar_t
                          ) / a_v[y, x]

            if x == n-1 and y == 1:
                a_v[y, x] = vv[y, x] * detar_x \
                            + 3 / 8 * detar_y * (uu[y, x+1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x]) \
                            + 3 / 2 * detar_y * detar_x / detar_t \
                            + 3 * detar_y / detar_x / Re \
                            + 2 * detar_x / detar_y / Re

                a_v[y, x+1] = 0
                a_v[y, x-1] = - 3 / 2 * detar_y / detar_x / Re \
                                - 3 / 8 * detar_y * (uu[y, x+1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x])
                a_v[y+1, x] = - detar_x / detar_y / Re
                a_v[y-1, x] = 0

                b = - 3 / 2 * vv[y, x] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[y, x] / alpha_v * vv[y, x]
                a_v[y, x] = -a_v[y, x] / alpha_v

                v[y, x] = (
                              a_v[y, x+1] * v[y, x+1]
                              + a_v[y, x-1] * v[y, x-1]
                              + a_v[y+1, x] * v[y+1, x]
                              + a_v[y-1, x] * v[y-1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y-1, x]) * detar_x
                          ) / a_v[y, x]

            if x == n-1 and y == n-1:
                a_v[y, x] = 3 / 2 * detar_y * detar_x \
                            + vv[y, x] * detar_x * detar_t \
                            + 3 / 8 * detar_y * detar_t * (uu[y, x+1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x]) \
                            + 3 * detar_y * detar_t / detar_x / Re \
                            + 2 * detar_x * detar_t / detar_y / Re

                a_v[y, x+1] = 0
                a_v[y, x-1] = - 3 / 2 * detar_y * detar_t / detar_x / Re \
                                - 3 / 8 * detar_y * detar_t * (
                                    uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x])
                a_v[y+1, x] = 0
                a_v[y-1, x] = - detar_x * detar_t / detar_y / Re \
                                - vv[y, x] * detar_x * detar_t

                b = - 3 / 2 * vv[y, x] * detar_x * detar_y

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[y, x] / alpha_v * vv[y, x]
                a_v[y, x] = -a_v[y, x] / alpha_v

                v[y, x] = (
                              a_v[y, x+1] * v[y, x+1]
                              + a_v[y, x-1] * v[y, x-1]
                              + a_v[y+1, x] * v[y+1, x]
                              + a_v[y-1, x] * v[y-1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y-1, x]) * detar_x * detar_t
                          ) / a_v[y, x]

            '四个边界'
            if x == 0 and y != 0 and y != 1 and y != n - 1 and y != n:
                a_v[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x]) / 4 * detar_y * detar_t \
                            + vv[y, x] * detar_x * detar_t \
                            + detar_x * detar_y \
                            + detar_y * detar_t / detar_x * 2 / Re \
                            + detar_x * detar_t / detar_y * 2 / Re

                a_v[y, x + 1] = - detar_y * detar_t / detar_x / Re
                a_v[y, x - 1] = 0
                a_v[y+1, x] = - detar_x * detar_t / detar_y / Re
                a_v[y-1, x] = - detar_x * detar_t / detar_y / Re \
                                - vv[y, x] * detar_x * detar_t

                b = - vv[y, x] * detar_x * detar_y

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[y, x] / alpha_v * vv[y, x]
                a_v[y, x] = -a_v[y, x] / alpha_v

                v[y, x] = (
                              a_v[y, x + 1] * v[y, x + 1]
                              + a_v[y, x - 1] * v[y, x - 1]
                              + a_v[y+1, x] * v[y+1, x]
                              + a_v[y-1, x] * v[y-1, x]
                              + b
                              + (p[y, x] - p[y-1, x]) * detar_x * detar_t
                          ) / a_v[y, x]

            if x == n - 1 and y!= 0 and y!= 1 and y!= n-1 and y!= n:
                a_v[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x]) / 4 * detar_y * detar_t \
                            + vv[y, x] * detar_x * detar_t \
                            + detar_x * detar_y \
                            + detar_y * detar_t / detar_x * 2 / Re \
                            + detar_x * detar_t / detar_y * 2 / Re

                a_v[y, x + 1] = 0
                a_v[y, x - 1] = - detar_y * detar_t / detar_x / Re \
                                - 4 * detar_y * detar_t * (uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x])
                a_v[y+1, x] = - detar_x * detar_t / detar_y / Re
                a_v[y-1, x] = - detar_x * detar_t / detar_y / Re \
                                - vv[y, x] * detar_x * detar_t

                b = - vv[y, x] * detar_x * detar_y

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[y, x] / alpha_v * vv[y, x]
                a_v[y, x] = -a_v[y, x] / alpha_v

                v[y, x] = (
                              a_v[y, x + 1] * v[y, x + 1]
                              + a_v[y, x - 1] * v[y, x - 1]
                              + a_v[y+1, x] * v[y+1, x]
                              + a_v[y-1, x] * v[y-1, x]
                              + b
                              + (p[y, x] - p[y-1, x]) * detar_x * detar_t
                          ) / a_v[y, x]

            if y == 1 and x != 0 and x != n - 1 and x != n:
                a_v[y, x] = vv[y, x] * detar_x \
                            + 3 / 8 * detar_y * (uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x]) \
                            + 3 / 2 * detar_y * detar_x / detar_t \
                            + 3 * detar_y / detar_x / Re \
                            + 2 * detar_x / detar_y / Re

                a_v[y, x + 1] = - 3 / 2 * detar_y / detar_x / Re
                a_v[y, x - 1] = - 3 / 2 * detar_y / detar_x / Re \
                                - 3 / 8 * detar_y * (uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x])
                a_v[y+1, x] = - detar_x / detar_y / Re
                a_v[y-1, x] = 0

                b = - 3 / 2 * vv[y, x] * detar_x * detar_y / detar_t

                '亚松弛设置'
                b = b + (1 - alpha_v) * a_v[y, x] / alpha_v * vv[y, x]
                a_v[y, x] = -a_v[y, x] / alpha_v

                v[y, x] = (
                              a_v[y, x + 1] * v[y, x + 1]
                              + a_v[y, x - 1] * v[y, x - 1]
                              + a_v[y+1, x] * v[y+1, x]
                              + a_v[y-1, x] * v[y-1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y-1, x]) * detar_x
                          ) / a_v[y, x]

            if y == n - 1 and x!=0 and x!=n-1 and x!=n:
                a_v[y, x] = 3 / 2 * detar_y * detar_x \
                            + vv[y, x] * detar_x * detar_t \
                            + 3 / 8 * detar_y * detar_t * (uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x]) \
                            + 3 * detar_y * detar_t / detar_x / Re \
                            + 2 * detar_x * detar_t / detar_y / Re

                a_v[y, x + 1] = - 3 / 2 * detar_y * detar_t / detar_x / Re
                a_v[y, x - 1] = - 3 / 2 * detar_y * detar_t / detar_x / Re \
                                - 3 / 8 * detar_y * detar_t * (uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x])
                a_v[y+1, x] = 0
                a_v[y-1, x] = - detar_x * detar_t / detar_y / Re \
                                - vv[y, x] * detar_x * detar_t

                b = - 3 / 2 * vv[y, x] * detar_x * detar_y

                '亚松弛设置'
                b = b + (1 - alpha_v)*a_v[y, x]/alpha_v*vv[y, x]
                a_v[y, x] = -a_v[y, x] / alpha_v

                v[y, x] = (
                              a_v[y, x + 1] * v[y, x + 1]
                              + a_v[y, x - 1] * v[y, x - 1]
                              + a_v[y+1, x] * v[y+1, x]
                              + a_v[y-1, x] * v[y-1, x]
                              + b
                              + 3 / 2 * (p[y, x] - p[y-1, x]) * detar_x * detar_t
                          ) / a_v[y, x]

            '内部节点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != 1 and y != n - 1 and y != n:
                a_v[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x]) / 4 * detar_y * detar_t \
                            + vv[y, x] * detar_x * detar_t \
                            + detar_x * detar_y \
                            + detar_y * detar_t / detar_x * 2 / Re \
                            + detar_x * detar_t / detar_y * 2 / Re

                a_v[y, x + 1] = - detar_y * detar_t / detar_x / Re
                a_v[y, x - 1] = - detar_y * detar_t / detar_x / Re \
                                - 4 * detar_y * detar_t * (uu[y, x + 1] + uu[y, x] + uu[y - 1, x + 1] + uu[y - 1, x])
                a_v[y+1, x] = - detar_x * detar_t / detar_y / Re
                a_v[y-1, x] = - detar_x * detar_t / detar_y / Re \
                                - vv[y, x] * detar_x * detar_t

                b = - vv[y, x] * detar_x * detar_y

                '亚松弛设置'
                b = b + (1 - alpha_v)*a_v[y, x]/alpha_v*vv[y, x]
                a_v[y, x] = -a_v[y, x] / alpha_v

                v[y, x] = (
                              a_v[y, x + 1] * v[y, x + 1]
                              + a_v[y, x - 1] * v[y, x - 1]
                              + a_v[y+1, x] * v[y+1, x]
                              + a_v[y-1, x] * v[y-1, x]
                              + b
                              + (p[y, x] - p[y-1, x]) * detar_x * detar_t
                          ) / a_v[y, x]
            "------------------------------------------压力修正方程计算---------------------------------------------"
    for y in range(n + 1):
        for x in range(n+1):
            '角点'
            if x == 0 and y == 0:
                a_u_add[y, x+1] = (vv[y, x+1] + vv[y, x] + vv[y+1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x+1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[y, x] = 1

                a_v_add[y+1, x] = (uu[y + 1, x + 1] + uu[y+1, x] + uu[y, x] + uu[y, x+1]) / 4 * detar_y * detar_t \
                                    + vv[y+1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re
                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y, x] = 1

                d[y, x+1] = detar_y * alpha_p / a_u_add[y, x+1]
                a_p[y, x+1] = d[y, x+1] * detar_y

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，aw为0,as=0'
                d[y, x-1] = 0
                a_p[y, x-1] = d[y, x-1] * detar_y

                d[y+1, x] = detar_x * alpha_p / a_v_add[y+1, x]
                a_p[y+1, x] = d[y+1, x] * detar_x

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，aw为0,as=0'
                d[y-1, x] = 0
                a_p[y-1, x] = d[y-1, x] * detar_x

                a_p[y, x] = a_p[y, x+1] + a_p[y, x-1] + a_p[y+1, x] + a_p[y-1, x]
                b_add = (u[y, x] - u[y, x+1]) * detar_y + (v[y, x] - v[y+1, x]) * detar_x

                p_add[y, x] = (
                                  a_p[y, x+1] * p_add[y, x+1] + a_p[y, x-1] * p_add[y, x-1] + a_p[y+1, x] *
                                  p_add[
                                      y+1, x] + a_p[y-1, x] * p_add[y-1, x] + b_add) / a_p[y, x]

                p[y, x] = p[y, x] + alpha_p * p_add[y, x]

            if x == n - 1 and y == 0:
                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[y, x+1] = 1

                a_u_add[y, x] = (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re


                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y, x] = 1

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，ae为0,as为0'
                d[y, x + 1] = 0
                a_p[y, x + 1] = d[y, x + 1] * detar_y

                d[y, x - 1] = detar_y * alpha_p/ a_u_add[y, x]
                a_p[y, x - 1] = d[y, x - 1] * detar_y

                d[y + 1, x] = detar_x * alpha_p / a_v_add[y + 1, x]
                a_p[y + 1, x] = d[y + 1, x] * detar_x

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，ae为0,as为0'
                d[y - 1, x] = 0
                a_p[y - 1, x] = d[y - 1, x] * detar_x

                a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
                b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

                p_add[y, x] = (
                                  a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1] + a_p[y + 1, x] *
                                  p_add[
                                      y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add) / a_p[y, x]

                p[y, x] = p[y, x] + alpha_p * p_add[y, x]

            if x == 0 and y == n - 1:
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re


                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[y, x] = 1

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y+1, x] = 1

                a_v_add[y, x] = (uu[y, x+1] + uu[y, x] + uu[y-1, x] + uu[y - 1, x+ 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[y, x + 1] = detar_y * alpha_p/ a_u_add[y, x+1]
                a_p[y, x + 1] = d[y, x + 1] * detar_y

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，aw为0，an为0'
                d[y, x - 1] = 0
                a_p[y, x - 1] = d[y, x - 1] * detar_y

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，aw为0，an为0'
                d[y + 1, x] = 0
                a_p[y + 1, x] = d[y + 1, x] * detar_x

                d[y-1, x] = detar_x * alpha_p/ a_v_add[y, x]
                a_p[y - 1, x] = d[y - 1, x] * detar_x

                a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
                b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

                p_add[y, x] = (
                                  a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1] + a_p[y + 1, x] *
                                  p_add[
                                      y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add) / a_p[y, x]

                p[y, x] = p[y, x] + alpha_p * p_add[y, x]

            if x == n - 1 and y == n - 1:
                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[y, x + 1] = 1

                a_u_add[y, x] = (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y+1, x] = 1

                a_v_add[y, x] = (uu[y, x+1] + uu[y, x] + uu[y-1, x] + uu[y - 1, x+ 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，an为0,ae为0'
                d[y, x+1] = 0
                a_p[y, x+1] = d[y, x+1] * detar_y

                d[y, x-1] = detar_y * alpha_p/ a_u_add[y, x]
                a_p[y, x-1] = d[y, x-1] * detar_y

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，an为0,ae为0'
                d[y+1, x] = 0
                a_p[y+1, x] = d[y+1, x] * detar_x

                d[y-1, x] = detar_x * alpha_p/ a_v_add[y, x]
                a_p[y-1, x] = d[y-1, x] * detar_x

                a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
                b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

                p_add[y, x] = (
                                  a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1] + a_p[y + 1, x] *
                                  p_add[
                                      y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add) / a_p[y, x]

                p[y, x] = p[y, x] + alpha_p * p_add[y, x]

            '边界节点'
            if x == 0 and y != 0 and y != n - 1 and y != n:
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左边界时，左边的速度u(i,j)为0，系数不使用置为１'
                a_u_add[y, x] = 1

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x] + uu[y - 1, x + 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[y, x + 1] = detar_y * alpha_p / a_u_add[y, x + 1]
                a_p[y, x + 1] = d[y, x + 1] * detar_y

                '左边界时，左边的速度u(i,j)为0，aw为0'
                d[y, x - 1] = 0
                a_p[y, x - 1] = d[y, x - 1] * detar_y

                d[y + 1, x] = detar_x * alpha_p / a_v_add[y + 1, x]
                a_p[y + 1, x] = d[y + 1, x] * detar_x

                d[y - 1, x] = detar_x * alpha_p / a_v_add[y, x]
                a_p[y - 1, x] = d[y - 1, x] * detar_x

                a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
                b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

                p_add[y, x] = (
                                  a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1] + a_p[y + 1, x] *
                                  p_add[
                                      y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add) / a_p[y, x]

                p[y, x] = p[y, x] + alpha_p * p_add[y, x]

            if x == n - 1 and y != 0 and y != n - 1 and y != n:
                '右边界时，右边的速度u(i,j)为0，系数不使用置为１'
                a_u_add[y, x + 1] = 1

                a_u_add[y, x] = (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x] + uu[y - 1, x + 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                '右边界时，右边的速度u(i,j)为0，ae系数为0'
                d[y, x + 1] = 0
                a_p[y, x + 1] = d[y, x + 1] * detar_y

                d[y, x - 1] = detar_y * alpha_p / a_u_add[y, x - 1]
                a_p[y, x - 1] = d[y, x - 1] * detar_y

                d[y + 1, x] = detar_x * alpha_p / a_v_add[y + 1, x]
                a_p[y + 1, x] = d[y + 1, x] * detar_x

                d[y - 1, x] = detar_x * alpha_p / a_v_add[y, x]
                a_p[y - 1, x] = d[y - 1, x] * detar_x

                a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
                b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

                p_add[y, x] = (
                                  a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1] + a_p[y + 1, x] *
                                  p_add[
                                      y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add) / a_p[y, x]

                p[y, x] = p[y, x] + alpha_p * p_add[y, x]

            if y == 0 and x != 0 and x != n - 1 and x != n:
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[y, x] = (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                '下边界时，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y, x] = 1

                d[y, x + 1] = detar_y * alpha_p / a_u_add[y, x + 1]
                a_p[y, x + 1] = d[y, x + 1] * detar_y

                d[y, x - 1] = detar_y * alpha_p / a_u_add[y, x - 1]
                a_p[y, x - 1] = d[y, x - 1] * detar_y

                d[y + 1, x] = detar_x * alpha_p / a_v_add[y + 1, x]
                a_p[y + 1, x] = d[y + 1, x] * detar_x

                '下边界时，下边的速度v(i,j)为0，as为0'
                d[y - 1, x] = 0
                a_p[y - 1, x] = d[y - 1, x] * detar_x

                a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
                b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

                p_add[y, x] = (
                                  a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1] + a_p[y + 1, x] *
                                  p_add[
                                      y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add) / a_p[y, x]

                p[y, x] = p[y, x] + alpha_p * p_add[y, x]

            if y == n - 1 and x != 0 and x != n - 1 and x != n:
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[y, x] = (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                '上边界时，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y+1, x] = 1

                a_v_add[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x] + uu[y - 1, x + 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[y, x + 1] = detar_y * alpha_p / a_u_add[y, x + 1]
                a_p[y, x + 1] = d[y, x + 1] * detar_y

                d[y, x - 1] = detar_y * alpha_p / a_u_add[y, x - 1]
                a_p[y, x - 1] = d[y, x - 1] * detar_y

                '上边界时，上边的速度v(i,j)为0，an系数为０'
                d[y + 1, x] = 0
                a_p[y + 1, x] = d[y + 1, x] * detar_x

                d[y - 1, x] = detar_x * alpha_p / a_v_add[y, x]
                a_p[y - 1, x] = d[y - 1, x] * detar_x

                a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
                b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

                p_add[y, x] = (
                                  a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1] + a_p[y + 1, x] *
                                  p_add[
                                      y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add) / a_p[y, x]

                p[y, x] = p[y, x] + alpha_p * p_add[y, x]

            '内部节点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != n - 1 and y != n:
                '补充方程系数'
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[y, x] = (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x] + uu[y - 1, x + 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[y, x + 1] = detar_y * alpha_p / a_u_add[y, x + 1]
                a_p[y, x + 1] = d[y, x + 1] * detar_y

                d[y, x - 1] = detar_y * alpha_p / a_u_add[y, x-1]
                a_p[y, x - 1] = d[y, x - 1] * detar_y

                d[y+1, x] = detar_x * alpha_p / a_v_add[y+1, x]
                a_p[y+1, x] = d[y+1, x] * detar_x

                d[y-1, x] = detar_x * alpha_p / a_v_add[y, x]
                a_p[y-1, x] = d[y-1, x] * detar_x

                a_p[y, x] = a_p[y, x + 1] + a_p[y, x - 1] + a_p[y + 1, x] + a_p[y - 1, x]
                b_add = (u[y, x] - u[y, x + 1]) * detar_y + (v[y, x] - v[y + 1, x]) * detar_x

                p_add[y, x] = (
                                a_p[y, x + 1] * p_add[y, x + 1] + a_p[y, x - 1] * p_add[y, x - 1]
                                + a_p[y + 1, x] *p_add[y + 1, x] + a_p[y - 1, x] * p_add[y - 1, x] + b_add
                              ) / a_p[y, x]

                p[y, x] = p[y, x] + alpha_p * p_add[y, x]
    "------------------------------------------速度修正计算---------------------------------------------"
    for y in range(n + 1):
        for x in range(n + 1):
            '角点'
            if x == 0 and y == 0:
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[y, x] = 1

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re
                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y, x] = 1

                d[y, x + 1] = detar_y * alpha_p / a_u_add[y, x + 1]

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，aw为0,as=0'
                d[y, x-1] = 0

                d[y+1, x] = detar_x * alpha_p/ a_v_add[y+1, x]

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，aw为0,as=0'
                d[y-1, x] = 0

                u[y, x + 1] = u[y, x + 1] + d[y, x + 1] * (p_add[y, x] - p_add[y, x + 1])
                v[y + 1, x] = v[y + 1, x] + d[y + 1, x] * (p_add[y, x] - p_add[y + 1, x])
                u[y, x - 1] = u[y, x - 1] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
                v[y - 1, x] = v[y - 1, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

            if x == n - 1 and y == 0:
                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[y, x + 1] = 1

                a_u_add[y, x] = (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y, x] = 1

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，ae为0,as为0'
                d[y, x+1] = 0

                d[y, x - 1] = detar_y * alpha_p / a_u_add[y, x]

                d[y + 1, x] = detar_x * alpha_p / a_v_add[y + 1, x]

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，ae为0,as为0'
                d[y - 1, x] = 0

                u[y, x + 1] = u[y, x + 1] + d[y, x + 1] * (p_add[y, x] - p_add[y, x + 1])
                v[y + 1, x] = v[y + 1, x] + d[y + 1, x] * (p_add[y, x] - p_add[y + 1, x])
                u[y, x - 1] = u[y, x - 1] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
                v[y - 1, x] = v[y - 1, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

            if x == 0 and y == n - 1:
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[y, x] = 1

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y + 1, x] = 1

                a_v_add[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x] + uu[y - 1, x + 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[y, x + 1] = detar_y * alpha_p / a_u_add[y, x + 1]

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，aw为0，an为0'
                d[y, x-1] = 0

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，aw为0，an为0'
                d[y+1, x] = 0

                d[y - 1, x] = detar_x * alpha_p / a_v_add[y, x]

                u[y, x + 1] = u[y, x + 1] + d[y, x + 1] * (p_add[y, x] - p_add[y, x + 1])
                v[y + 1, x] = v[y + 1, x] + d[y + 1, x] * (p_add[y, x] - p_add[y + 1, x])
                u[y, x - 1] = u[y, x - 1] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
                v[y - 1, x] = v[y - 1, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

            if x == n - 1 and y == n - 1:
                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[y, x + 1] = 1

                a_u_add[y, x] = (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y + 1, x] = 1

                a_v_add[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x] + uu[y - 1, x + 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，an为0,ae为0'
                d[y, x + 1] = 0

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，an为0,ae为0'
                d[y + 1, x] = 0

                d[y, x - 1] = detar_y * alpha_p / a_u_add[y, x]

                d[y - 1, x] = detar_x * alpha_p / a_v_add[y, x]

                u[y, x + 1] = u[y, x + 1] + d[y, x + 1] * (p_add[y, x] - p_add[y, x + 1])
                v[y + 1, x] = v[y + 1, x] + d[y + 1, x] * (p_add[y, x] - p_add[y + 1, x])
                u[y, x - 1] = u[y, x - 1] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
                v[y - 1, x] = v[y - 1, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

            '边界节点'
            if x == 0 and y!= 0 and y!= n-1 and y!= n:
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左边界时，左边的速度u(i,j)为0，系数不使用置为１'
                a_u_add[y, x] = 1

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x] + uu[y - 1, x + 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[y, x + 1] = detar_y * alpha_p / a_u_add[y, x + 1]

                d[y + 1, x] = detar_x * alpha_p / a_v_add[y + 1, x]

                '左边界时，左边的速度u(i,j)为0，aw为0'
                d[y, x - 1] = 0

                d[y - 1, x] = detar_x * alpha_p / a_v_add[y, x]

                u[y, x + 1] = u[y, x + 1] + d[y, x + 1] * (p_add[y, x] - p_add[y, x + 1])
                v[y + 1, x] = v[y + 1, x] + d[y + 1, x] * (p_add[y, x] - p_add[y + 1, x])
                u[y, x - 1] = u[y, x - 1] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
                v[y - 1, x] = v[y - 1, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

            if x == n-1 and y != 0 and y != n - 1 and y != n:
                '右边界时，右边的速度u(i,j)为0，系数不使用置为１'
                a_u_add[y, x + 1] = 1

                a_u_add[y, x] = (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x] + uu[y - 1, x + 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                '右边界时，右边的速度u(i,j)为0，ae系数为0'
                d[y, x + 1] = 0

                d[y + 1, x] = detar_x * alpha_p / a_v_add[y + 1, x]

                d[y, x - 1] = detar_y * alpha_p / a_u_add[y, x]

                d[y - 1, x] = detar_x * alpha_p / a_v_add[y, x]

                u[y, x + 1] = u[y, x + 1] + d[y, x + 1] * (p_add[y, x] - p_add[y, x + 1])
                v[y + 1, x] = v[y + 1, x] + d[y + 1, x] * (p_add[y, x] - p_add[y + 1, x])
                u[y, x - 1] = u[y, x - 1] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
                v[y - 1, x] = v[y - 1, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

            if y == 0 and x != 0 and x != n - 1 and x != n:
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[y, x] = (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[y + 1, x] = (uu[y + 1, x + 1] + uu[y + 1, x] + uu[y, x] + uu[y, x + 1]) / 4 * detar_y * detar_t \
                                    + vv[y + 1, x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                '下边界时，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y, x] = 1

                d[y, x + 1] = detar_y * alpha_p / a_u_add[y, x + 1]

                d[y + 1, x] = detar_x * alpha_p / a_v_add[y + 1, x]

                d[y, x - 1] = detar_y * alpha_p / a_u_add[y, x]

                '下边界时，下边的速度v(i,j)为0，as为0'
                d[y - 1, x] = 0

                u[y, x + 1] = u[y, x + 1] + d[y, x + 1] * (p_add[y, x] - p_add[y, x + 1])
                v[y + 1, x] = v[y + 1, x] + d[y + 1, x] * (p_add[y, x] - p_add[y + 1, x])
                u[y, x - 1] = u[y, x - 1] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
                v[y - 1, x] = v[y - 1, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

            if y == n-1 and x != 0 and x != n - 1 and x != n:
                a_u_add[y, x + 1] = (vv[y, x + 1] + vv[y, x] + vv[y + 1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x + 1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[y, x] = (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                '上边界时，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[y + 1, x] = 1

                a_v_add[y, x] = (uu[y, x + 1] + uu[y, x] + uu[y - 1, x] + uu[y - 1, x + 1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[y, x + 1] = detar_y * alpha_p / a_u_add[y, x + 1]

                '上边界时，上边的速度v(i,j)为0，an系数为０'
                d[y + 1, x] = 0

                d[y, x - 1] = detar_y * alpha_p / a_u_add[y, x]

                d[y - 1, x] = detar_x * alpha_p / a_v_add[y, x]

                u[y, x + 1] = u[y, x + 1] + d[y, x + 1] * (p_add[y, x] - p_add[y, x + 1])
                v[y + 1, x] = v[y + 1, x] + d[y + 1, x] * (p_add[y, x] - p_add[y + 1, x])
                u[y, x - 1] = u[y, x - 1] + d[y, x - 1] * (p_add[y, x - 1] - p_add[y, x])
                v[y - 1, x] = v[y - 1, x] + d[y - 1, x] * (p_add[y - 1, x] - p_add[y, x])

            '内部节点'
            if x != 0 and x != n - 1 and x != n and y != 0 and y != n - 1 and y != n:
                a_u_add[y, x+1] = (vv[y, x+1] + vv[y, x] + vv[y+1, x] + vv[y + 1, x + 1]) / 4 * detar_x \
                                    + uu[y, x+1] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[y, x] = (vv[y, x] + vv[y, x-1] + vv[y+1, x-1] + vv[y+1, x]) / 4 * detar_x \
                                + uu[y, x] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[y+1,x] = (uu[y+1,x+1] + uu[y+1,x] + uu[y,x] + uu[y,x+1]) / 4 * detar_y * detar_t \
                                    + vv[y+1,x] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[y, x] = (uu[y, x+1] + uu[y, x] + uu[y-1, x] + uu[y-1, x+1]) / 4 * detar_y * detar_t \
                                + vv[y, x] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[y, x+1] = detar_y* alpha_p/a_u_add[y, x+1]

                d[y+1, x]= detar_x* alpha_p/a_v_add[y+1, x]

                d[y, x-1]= detar_y* alpha_p/a_u_add[y, x]

                d[y-1, x]= detar_x* alpha_p/a_v_add[y, x]

                u[y, x+1] = u[y, x+1] + d[y, x+1]*(p_add[y, x] - p_add[y, x+1])
                v[y+1, x] = v[y+1, x] + d[y+1, x]*(p_add[y, x] - p_add[y+1, x])
                u[y, x-1] = u[y, x-1] + d[y, x-1]*(p_add[y, x-1] - p_add[y, x])
                v[y-1, x] = v[y-1, x] + d[y-1, x]*(p_add[y-1, x] - p_add[y, x])

    '将本轮迭代的u,v保存'
    # vv = v
    # uu = u
    print(u)
    # print(v)