import numpy as np

"---------------------------------------------Problem-------------------------------------------------------------------"

'假设0＜=x,y<=1的方腔内充满粘性不可压缩流体，左右，下壁固定，上壁以u=-16x^2(1-x^2)运动，试求Re=100,200,400时的定常解'

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
alpha_u = 0.7
alpha_v = 0.7
alpha_p = 0.3
"------------------------------------------Iterative_Solution----------------------------------------------------------"

u[5] = 1
print(u)

for k in range(60):
    for y in range(n+1):
        for x in range(n+1):
            "------------------------------------------u控制方程计算---------------------------------------------"
            '四个脚点'
            if i == 0 and j == 1:
                a_u[y, x] = 3 / 2 * detar_x * detar_y / detar_t \
                            + uu[y, x] * detar_y \
                            + 3 * detar_x / Re / detar_y \
                            + 3 / 8 * detar_x * (vv[y, x] + vv[y, x - 1] + vv[y + 1, x - 1] + vv[y + 1, x]) \
                            + 2 / Re * detar_y / detar_x

                a_u[i + 1, j] = - detar_y / detar_x / Re
                a_u[i - 1, j] = 0
                a_u[i, j + 1] = - 3 / 2 * detar_x / detar_y / Re
                a_u[i, j - 1] = 0

                b = - 3 / 2 * uu[i, j] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[i, j] / alpha_u * uu[i, j]
                a_u[i, j] = -a_u[i, j] / alpha_u

                u[i, j] = (
                              a_u[i + 1, j] * u[i + 1, j]
                              + a_u[i - 1, j] * u[i - 1, j]
                              + a_u[i, j + 1] * u[i, j + 1]
                              + a_u[i, j - 1] * u[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i - 1, j]) * detar_y
                          ) / a_u[i, j]

            if i == 1 and j == n-1:
                a_u[i, j] = 3 / 2 * detar_x * detar_y / detar_t \
                            + uu[i, j] * detar_y \
                            + 3 * detar_x / Re / detar_y \
                            + 3 / 8 * detar_x * (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) \
                            + 2 / Re * detar_y / detar_x

                a_u[i + 1, j] = - detar_y / detar_x / Re
                a_u[i - 1, j] = 0
                a_u[i, j + 1] = 0
                a_u[i, j - 1] = - 3 / 8 * detar_x * (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) \
                                - 3 / 2 / Re * detar_x / detar_y

                b = - 3 / 2 * uu[i, j] * detar_x * detar_y / detar_t - 3/2*detar_x/Re/detar_y

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[i, j] / alpha_u * uu[i, j]
                a_u[i, j] = -a_u[i, j] / alpha_u

                u[i, j] = (
                              a_u[i + 1, j] * u[i + 1, j]
                              + a_u[i - 1, j] * u[i - 1, j]
                              + a_u[i, j + 1] * u[i, j + 1]
                              + a_u[i, j - 1] * u[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i - 1, j]) * detar_y
                          ) / a_u[i, j]

            if i == n-1 and j == 0:
                a_u[i, j] = 3 / 2 * detar_x * detar_y / detar_t \
                            + 3 / 8 * detar_x * (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 3 / 2 * detar_x / detar_y

                a_u[i + 1, j] = 0
                a_u[i - 1, j] = -uu[i, j] * detar_y - detar_y / Re / detar_x
                a_u[i, j + 1] = -3 / 2 * detar_x / detar_y / Re
                a_u[i, j - 1] = 0

                b = - 3 / 2 * uu[i, j] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[i, j] / alpha_u * uu[i, j]
                a_u[i, j] = -a_u[i, j] / alpha_u

                u[i, j] = (
                              a_u[i + 1, j] * u[i + 1, j]
                              + a_u[i - 1, j] * u[i - 1, j]
                              + a_u[i, j + 1] * u[i, j + 1]
                              + a_u[i, j - 1] * u[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i - 1, j]) * detar_y
                          ) / a_u[i, j]

            if i == n-1 and j == n-1:
                a_u[i, j] = 3 / 2 * detar_x * detar_y / detar_t \
                            + 3 / 8 * detar_x * (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 3 / 2 * detar_x / detar_y

                a_u[i + 1, j] = 0
                a_u[i - 1, j] = -uu[i, j] * detar_y - detar_y / Re / detar_x
                a_u[i, j + 1] = 0
                a_u[i, j - 1] = -3 / 8 * detar_x * (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) \
                                - 3 / 2 * detar_x / Re / detar_y

                b = - 3 / 2 * uu[i, j] * detar_x * detar_y / detar_t - 3 / 2 * detar_x / detar_y / Re

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[i, j] / alpha_u * uu[i, j]
                a_u[i, j] = -a_u[i, j] / alpha_u

                u[i, j] = (
                              a_u[i + 1, j] * u[i + 1, j]
                              + a_u[i - 1, j] * u[i - 1, j]
                              + a_u[i, j + 1] * u[i, j + 1]
                              + a_u[i, j - 1] * u[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i - 1, j]) * detar_y
                        ) / a_u[i, j]

            '四个边界'
            if i == 1 and j != n-1 and j !=0 and j!=n:
                a_u[i, j] = 3 / 2 * detar_x * detar_y / detar_t \
                            + uu[i, j] * detar_y \
                            + 3 * detar_x / Re / detar_y \
                            + 3 / 8 * detar_x * (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) \
                            + 2 / Re * detar_y / detar_x

                a_u[i + 1, j] = - detar_y / detar_x / Re
                a_u[i - 1, j] = 0
                a_u[i, j + 1] = - 3 / 2 * detar_x / detar_y / Re
                a_u[i, j - 1] = - 3 / 8 * detar_x * (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) \
                                - 3 / 2 / Re * detar_x / detar_y

                b = - 3 / 2 * uu[i, j] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[i, j] / alpha_u * uu[i, j]
                a_u[i, j] = -a_u[i, j] / alpha_u

                u[i, j] = (
                              a_u[i + 1, j] * u[i + 1, j]
                              + a_u[i - 1, j] * u[i - 1, j]
                              + a_u[i, j + 1] * u[i, j + 1]
                              + a_u[i, j - 1] * u[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i - 1, j]) * detar_y
                          ) / a_u[i, j]

            if i == n-1 and j!= n-1 and j != 0 and j!=n:
                a_u[i, j] = 3 / 2 * detar_x * detar_y / detar_t \
                            + 3 / 8 * detar_x * (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) \
                            + 2 / Re * detar_y / detar_x \
                            + 2 / Re * 3 / 2 * detar_x / detar_y

                a_u[i + 1, j] = 0
                a_u[i - 1, j] = -uu[i, j] * detar_y - detar_y / Re / detar_x
                a_u[i, j + 1] = -3 / 2 * detar_x / detar_y / Re
                a_u[i, j - 1] = -3 / 8 * detar_x * (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) \
                                - 3 / 2 * detar_x / Re / detar_y

                b = - 3 / 2 * uu[i, j] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[i, j] / alpha_u * uu[i, j]
                a_u[i, j] = -a_u[i, j] / alpha_u

                u[i, j] = (
                              a_u[i + 1, j] * u[i + 1, j]
                              + a_u[i - 1, j] * u[i - 1, j]
                              + a_u[i, j + 1] * u[i, j + 1]
                              + a_u[i, j - 1] * u[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i - 1, j]) * detar_y
                          ) / a_u[i, j]

            if j == 0 and i!=0 and i!= 1 and i != n-1 and i!=n:
                a_u[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                            + uu[i, j] * detar_y \
                            + detar_x * detar_y / detar_t \
                            + detar_y / detar_x * 2 / Re \
                            + detar_x / detar_y * 2 / Re

                a_u[i + 1, j] = - detar_y / detar_x / Re
                a_u[i - 1, j] = - uu[i, j] * detar_y \
                                - detar_y / detar_x / Re
                a_u[i, j + 1] = - detar_x / detar_y / Re
                a_u[i, j - 1] = 0

                b = - uu[i, j] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[i, j] / alpha_u * uu[i, j]
                a_u[i, j] = -a_u[i, j] / alpha_u

                u[i, j] = (
                              a_u[i + 1, j] * u[i + 1, j]
                              + a_u[i - 1, j] * u[i - 1, j]
                              + a_u[i, j + 1] * u[i, j + 1]
                              + a_u[i, j - 1] * u[i, j - 1]
                              + b
                              + (p[i, j] - p[i - 1, j]) * detar_y
                          ) / a_u[i, j]

            if j == n - 1 and i!=0 and i!= 1 and i !=n-1 and i!= n:
                a_u[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                            + uu[i, j] * detar_y \
                            + detar_x * detar_y / detar_t \
                            + detar_y / detar_x * 2 / Re \
                            + detar_x / detar_y * 2 / Re

                a_u[i + 1, j] = - detar_y / detar_x / Re
                a_u[i - 1, j] = - uu[i, j] * detar_y \
                                - detar_y / detar_x / Re
                a_u[i, j + 1] = 0
                a_u[i, j - 1] = - (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                - detar_x / Re / detar_y

                b = - uu[i, j] * detar_x * detar_y / detar_t - 3 / 2 * detar_x / detar_y / Re

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[i, j] / alpha_u * uu[i, j]
                a_u[i, j] = -a_u[i, j] / alpha_u

                u[i, j] = (
                              a_u[i + 1, j] * u[i + 1, j]
                              + a_u[i - 1, j] * u[i - 1, j]
                              + a_u[i, j + 1] * u[i, j + 1]
                              + a_u[i, j - 1] * u[i, j - 1]
                              + b
                              + (p[i, j] - p[i - 1, j]) * detar_y
                          ) / a_u[i, j]

            '内部节点'
            if i != 0 and i!= 1 and i != n - 1 and i!= n and j != 0 and j != n - 1 and j!= n:
                a_u[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                            + uu[i, j] * detar_y \
                            + detar_x * detar_y / detar_t \
                            + detar_y / detar_x * 2 / Re\
                            + detar_x / detar_y * 2 / Re

                a_u[i + 1, j] = - detar_y / detar_x / Re
                a_u[i - 1, j] = - uu[i, j] * detar_y \
                                - detar_y / detar_x / Re
                a_u[i, j + 1] = - detar_x / detar_y / Re
                a_u[i, j - 1] = - (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                - detar_x / Re / detar_y

                b = - uu[i, j] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_u) * a_u[i, j] / alpha_u * uu[i, j]
                a_u[i, j] = -a_u[i, j] / alpha_u

                u[i, j] = (
                              a_u[i + 1, j] * u[i + 1, j]
                              + a_u[i - 1, j] * u[i - 1, j]
                              + a_u[i, j + 1] * u[i, j + 1]
                              + a_u[i, j - 1] * u[i, j - 1]
                              + b
                              + (p[i, j] - p[i - 1, j]) * detar_y
                          ) / a_u[i, j]


            "------------------------------------------v控制方程计算---------------------------------------------"
            '四个脚点'
            if i == 0 and j == 1:
                a_v[i, j] = vv[i, j] * detar_x \
                            + 3 / 8 * detar_y * (uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1]) \
                            + 3 / 2 * detar_y * detar_x / detar_t \
                            + 3 * detar_y / detar_x / Re \
                            + 2 * detar_x / detar_y / Re

                a_v[i + 1, j] = - 3 / 2 * detar_y / detar_x / Re
                a_v[i - 1, j] = 0
                a_v[i, j + 1] = - detar_x / detar_y / Re
                a_v[i, j - 1] = 0

                b = - 3 / 2 * vv[i, j] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[i, j] / alpha_v * vv[i, j]
                a_v[i, j] = -a_v[i, j] / alpha_v

                v[i, j] = (
                              a_v[i + 1, j] * v[i + 1, j]
                              + a_v[i - 1, j] * v[i - 1, j]
                              + a_v[i, j + 1] * v[i, j + 1]
                              + a_v[i, j - 1] * v[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i, j - 1]) * detar_x
                          ) / a_v[i, j]

            if i == 0 and j == n-1:
                a_v[i, j] = 3 / 2 * detar_y * detar_x \
                            + vv[i, j] * detar_x * detar_t \
                            + 3 / 8 * detar_y * detar_t * (uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1]) \
                            + 3 * detar_y * detar_t / detar_x / Re \
                            + 2 * detar_x * detar_t / detar_y / Re

                a_v[i + 1, j] = - 3 / 2 * detar_y * detar_t / detar_x / Re
                a_v[i - 1, j] = 0
                a_v[i, j + 1] = 0
                a_v[i, j - 1] = - detar_x * detar_t / detar_y / Re \
                                - vv[i, j] * detar_x * detar_t

                b = - 3 / 2 * vv[i, j] * detar_x * detar_y

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[i, j] / alpha_v * vv[i, j]
                a_v[i, j] = -a_v[i, j] / alpha_v

                v[i, j] = (
                              a_v[i + 1, j] * v[i + 1, j]
                              + a_v[i - 1, j] * v[i - 1, j]
                              + a_v[i, j + 1] * v[i, j + 1]
                              + a_v[i, j - 1] * v[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i, j - 1]) * detar_x * detar_t
                          ) / a_v[i, j]

            if i == n-1 and j == 1:
                a_v[i, j] = vv[i, j] * detar_x \
                            + 3 / 8 * detar_y * (uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1]) \
                            + 3 / 2 * detar_y * detar_x / detar_t \
                            + 3 * detar_y / detar_x / Re \
                            + 2 * detar_x / detar_y / Re

                a_v[i + 1, j] = 0
                a_v[i - 1, j] = - 3 / 2 * detar_y / detar_x / Re \
                                - 3 / 8 * detar_y * (uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1])
                a_v[i, j + 1] = - detar_x / detar_y / Re
                a_v[i, j - 1] = 0

                b = - 3 / 2 * vv[i, j] * detar_x * detar_y / detar_t

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[i, j] / alpha_v * vv[i, j]
                a_v[i, j] = -a_v[i, j] / alpha_v

                v[i, j] = (
                              a_v[i + 1, j] * v[i + 1, j]
                              + a_v[i - 1, j] * v[i - 1, j]
                              + a_v[i, j + 1] * v[i, j + 1]
                              + a_v[i, j - 1] * v[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i, j - 1]) * detar_x
                          ) / a_v[i, j]

            if i == n-1 and j == n-1:
                a_v[i, j] = 3 / 2 * detar_y * detar_x \
                            + vv[i, j] * detar_x * detar_t \
                            + 3 / 8 * detar_y * detar_t * (uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1]) \
                            + 3 * detar_y * detar_t / detar_x / Re \
                            + 2 * detar_x * detar_t / detar_y / Re

                a_v[i + 1, j] = 0
                a_v[i - 1, j] = - 3 / 2 * detar_y * detar_t / detar_x / Re \
                                - 3 / 8 * detar_y * detar_t * (
                                uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1])
                a_v[i, j + 1] = 0
                a_v[i, j - 1] = - detar_x * detar_t / detar_y / Re \
                                - vv[i, j] * detar_x * detar_t

                b = - 3 / 2 * vv[i, j] * detar_x * detar_y

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[i, j] / alpha_v * vv[i, j]
                a_v[i, j] = -a_v[i, j] / alpha_v

                v[i, j] = (
                              a_v[i + 1, j] * v[i + 1, j]
                              + a_v[i - 1, j] * v[i - 1, j]
                              + a_v[i, j + 1] * v[i, j + 1]
                              + a_v[i, j - 1] * v[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i, j - 1]) * detar_x * detar_t
                          ) / a_v[i, j]

            # '四个边界'
            if i == 0 and j != 0 and j !=1 and j!= n-1 and j!=n :
                a_v[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                            + vv[i, j] * detar_x * detar_t \
                            + detar_x * detar_y \
                            + detar_y * detar_t / detar_x * 2 / Re \
                            + detar_x * detar_t / detar_y * 2 / Re

                a_v[i + 1, j] = - detar_y * detar_t / detar_x / Re
                a_v[i - 1, j] = 0
                a_v[i, j + 1] = - detar_x * detar_t / detar_y / Re
                a_v[i, j - 1] = - detar_x * detar_t / detar_y / Re \
                                - vv[i, j] * detar_x * detar_t

                b = - vv[i, j] * detar_x * detar_y

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[i, j] / alpha_v * vv[i, j]
                a_v[i, j] = -a_v[i, j] / alpha_v

                v[i, j] = (
                              a_v[i + 1, j] * v[i + 1, j]
                              + a_v[i - 1, j] * v[i - 1, j]
                              + a_v[i, j + 1] * v[i, j + 1]
                              + a_v[i, j - 1] * v[i, j - 1]
                              + b
                              + (p[i, j] - p[i, j - 1]) * detar_x * detar_t
                          ) / a_v[i, j]

            if i == n - 1 and j!= 0 and j!= 1 and j!= n-1 and j!= n:
                a_v[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                            + vv[i, j] * detar_x * detar_t \
                            + detar_x * detar_y \
                            + detar_y * detar_t / detar_x * 2 / Re \
                            + detar_x * detar_t / detar_y * 2 / Re

                a_v[i + 1, j] = 0
                a_v[i - 1, j] = - detar_y * detar_t / detar_x / Re \
                                - 4 * detar_y * detar_t * (uu[i + 1, j] + uu[i, j] + uu[i - 1, j - 1] + uu[i - 1, j])
                a_v[i, j + 1] = - detar_x * detar_t / detar_y / Re
                a_v[i, j - 1] = - detar_x * detar_t / detar_y / Re \
                                - vv[i, j] * detar_x * detar_t

                b = - vv[i, j] * detar_x * detar_y

                '亚松弛设置,注意两者的顺序'
                b = b + (1 - alpha_v) * a_v[i, j] / alpha_v * vv[i, j]
                a_v[i, j] = -a_v[i, j] / alpha_v

                v[i, j] = (
                              a_v[i + 1, j] * v[i + 1, j]
                              + a_v[i - 1, j] * v[i - 1, j]
                              + a_v[i, j + 1] * v[i, j + 1]
                              + a_v[i, j - 1] * v[i, j - 1]
                              + b
                              + (p[i, j] - p[i, j - 1]) * detar_x * detar_t
                          ) / a_v[i, j]

            if j == 1 and i!=0 and i!=n-1 and i!=n:
                a_v[i, j] = vv[i, j] * detar_x \
                            + 3 / 8 * detar_y * (uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1]) \
                            + 3 / 2 * detar_y * detar_x / detar_t \
                            + 3 * detar_y / detar_x / Re \
                            + 2 * detar_x / detar_y / Re

                a_v[i + 1, j] = - 3 / 2 * detar_y / detar_x / Re
                a_v[i - 1, j] = - 3 / 2 * detar_y / detar_x / Re \
                                - 3 / 8 * detar_y * (uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1])
                a_v[i, j + 1] = - detar_x / detar_y / Re
                a_v[i, j - 1] = 0

                b = - 3 / 2 * vv[i, j] * detar_x * detar_y / detar_t

                '亚松弛设置'
                b = b + (1 - alpha_v) * a_v[i, j] / alpha_v * vv[i, j]
                a_v[i, j] = -a_v[i, j] / alpha_v

                v[i, j] = (
                              a_v[i + 1, j] * v[i + 1, j]
                              + a_v[i - 1, j] * v[i - 1, j]
                              + a_v[i, j + 1] * v[i, j + 1]
                              + a_v[i, j - 1] * v[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i, j - 1]) * detar_x
                          ) / a_v[i, j]

            if j == n - 1 and i!=0 and i!=n-1 and i!=n:
                a_v[i, j] = 3 / 2 * detar_y * detar_x \
                            + vv[i, j] * detar_x * detar_t \
                            + 3 / 8 * detar_y * detar_t * (uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1]) \
                            + 3 * detar_y * detar_t / detar_x / Re \
                            + 2 * detar_x * detar_t / detar_y / Re

                a_v[i + 1, j] = - 3 / 2 * detar_y * detar_t / detar_x / Re
                a_v[i - 1, j] = - 3 / 2 * detar_y * detar_t / detar_x / Re \
                                - 3 / 8 * detar_y * detar_t * (uu[i + 1, j] + uu[i, j] + uu[i + 1, j - 1] + uu[i, j - 1])
                a_v[i, j + 1] = 0
                a_v[i, j - 1] = - detar_x * detar_t / detar_y / Re \
                                - vv[i, j] * detar_x * detar_t

                b = - 3 / 2 * vv[i, j] * detar_x * detar_y

                '亚松弛设置'
                b = b + (1 - alpha_v)*a_v[i,j]/alpha_v*vv[i,j]
                a_v[i, j] = -a_v[i, j] / alpha_v

                v[i, j] = (
                              a_v[i + 1, j] * v[i + 1, j]
                              + a_v[i - 1, j] * v[i - 1, j]
                              + a_v[i, j + 1] * v[i, j + 1]
                              + a_v[i, j - 1] * v[i, j - 1]
                              + b
                              + 3 / 2 * (p[i, j] - p[i, j - 1]) * detar_x * detar_t
                          ) / a_v[i, j]

            '内部节点'
            if i != 0 and i != n - 1 and i != n and j != 0 and j != 1 and j != n - 1 and j != n:
                a_v[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                            + vv[i, j] * detar_x * detar_t \
                            + detar_x * detar_y \
                            + detar_y * detar_t / detar_x * 2 / Re \
                            + detar_x * detar_t / detar_y * 2 / Re

                a_v[i + 1, j] = - detar_y * detar_t / detar_x / Re
                a_v[i - 1, j] = - detar_y * detar_t / detar_x / Re \
                                - 4 * detar_y * detar_t * (uu[i + 1, j] + uu[i, j] + uu[i - 1, j - 1] + uu[i - 1, j])
                a_v[i, j + 1] = - detar_x * detar_t / detar_y / Re
                a_v[i, j - 1] = - detar_x * detar_t / detar_y / Re \
                                - vv[i, j] * detar_x * detar_t

                b = - vv[i, j] * detar_x * detar_y

                '亚松弛设置'
                b = b + (1 - alpha_v)*a_v[i,j]/alpha_v*vv[i,j]
                a_v[i, j] = -a_v[i, j] / alpha_v

                v[i, j] = (
                              a_v[i + 1, j] * v[i + 1, j]
                              + a_v[i - 1, j] * v[i - 1, j]
                              + a_v[i, j + 1] * v[i, j + 1]
                              + a_v[i, j - 1] * v[i, j - 1]
                              + b
                              + (p[i, j] - p[i, j - 1]) * detar_x * detar_t
                          ) / a_v[i, j]

            "------------------------------------------压力修正方程计算---------------------------------------------"
    for i in range(n + 1):
        for j in range(n+1):
            '角点'
            if i == 0 and j == 0:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[i, j] = 1

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re
                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j] = 1

                d[i + 1, j] = detar_y * alpha_p / a_u_add[i + 1, j]
                a_p[i + 1, j] = d[i + 1, j] * detar_y

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，aw为0,as=0'
                d[i - 1, j] = 0
                a_p[i - 1, j] = d[i - 1, j] * detar_y

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]
                a_p[i, j + 1] = d[i, j + 1] * detar_x

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，aw为0,as=0'
                d[i, j - 1] = 0
                a_p[i, j - 1] = d[i, j - 1] * detar_x


                a_p[i, j] = a_p[i + 1, j] + a_p[i - 1, j] + a_p[i, j + 1] + a_p[i, j - 1]
                b_add = (u[i, j] - u[i + 1, j]) * detar_y + (v[i, j] - v[i, j + 1]) * detar_x

                p_add[i, j] = (
                                  a_p[i + 1, j] * p_add[i + 1, j] + a_p[i - 1, j] * p_add[i - 1, j] + a_p[i, j + 1] *
                                  p_add[
                                      i, j + 1] + a_p[i, j - 1] * p_add[i, j - 1] + b_add) / a_p[i, j]

                p[i, j] = p[i, j] + alpha_p * p_add[i, j]

            if i == n - 1 and j == 0:
                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[i + 1, j] = 1

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j] = 1

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，ae为0,as为0'
                d[i + 1, j] = 0
                a_p[i + 1, j] = d[i + 1, j] * detar_y

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]
                a_p[i - 1, j] = d[i - 1, j] * detar_y

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]
                a_p[i, j + 1] = d[i, j + 1] * detar_x

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，ae为0,as为0'
                d[i, j - 1] = 0
                a_p[i, j - 1] = d[i, j - 1] * detar_x

                a_p[i, j] = a_p[i + 1, j] + a_p[i - 1, j] + a_p[i, j + 1] + a_p[i, j - 1]
                b_add = (u[i, j] - u[i + 1, j]) * detar_y + (v[i, j] - v[i, j + 1]) * detar_x

                p_add[i, j] = (
                                  a_p[i + 1, j] * p_add[i + 1, j] + a_p[i - 1, j] * p_add[i - 1, j] + a_p[i, j + 1] *
                                  p_add[
                                      i, j + 1] + a_p[i, j - 1] * p_add[i, j - 1] + b_add) / a_p[i, j]

                p[i, j] = p[i, j] + alpha_p * p_add[i, j]

            if i == 0 and j == n - 1:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[i, j] = 1

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j + 1] = 1

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[i + 1, j] = detar_y * alpha_p/ a_u_add[i + 1, j]
                a_p[i + 1, j] = d[i + 1, j] * detar_y

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，aw为0，an为0'
                d[i - 1, j] = 0
                a_p[i - 1, j] = d[i - 1, j] * detar_y

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，aw为0，an为0'
                d[i, j + 1] = 0
                a_p[i, j + 1] = d[i, j + 1] * detar_x

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]
                a_p[i, j - 1] = d[i, j - 1] * detar_x

                a_p[i, j] = a_p[i + 1, j] + a_p[i - 1, j] + a_p[i, j + 1] + a_p[i, j - 1]
                b_add = (u[i, j] - u[i + 1, j]) * detar_y + (v[i, j] - v[i, j + 1]) * detar_x

                p_add[i, j] = (
                                  a_p[i + 1, j] * p_add[i + 1, j] + a_p[i - 1, j] * p_add[i - 1, j] + a_p[i, j + 1] *
                                  p_add[i, j + 1] + a_p[i, j - 1] * p_add[i, j - 1] + b_add) / a_p[i, j]

                p[i, j] = p[i, j] + alpha_p * p_add[i, j]

            if i == n - 1 and j == n - 1:
                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[i + 1, j] = 1

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j + 1] = 1

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，an为0,ae为0'
                d[i + 1, j] = 0
                a_p[i + 1, j] = d[i + 1, j] * detar_y

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]
                a_p[i - 1, j] = d[i - 1, j] * detar_y

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，an为0,ae为0'
                d[i, j + 1] = 0
                a_p[i, j + 1] = d[i, j + 1] * detar_x

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]
                a_p[i, j - 1] = d[i, j - 1] * detar_x

                a_p[i, j] = a_p[i + 1, j] + a_p[i - 1, j] + a_p[i, j + 1] + a_p[i, j - 1]
                b_add = (u[i, j] - u[i + 1, j]) * detar_y + (v[i, j] - v[i, j + 1]) * detar_x

                p_add[i, j] = (
                                  a_p[i + 1, j] * p_add[i + 1, j] + a_p[i - 1, j] * p_add[i - 1, j] + a_p[i, j + 1] *
                                  p_add[i, j + 1] + a_p[i, j - 1] * p_add[i, j - 1] + b_add) / a_p[i, j]

                p[i, j] = p[i, j] + alpha_p * p_add[i, j]

            '边界节点'
            if i == 0 and j!= 0 and j!= n-1 and j!= n:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左边界时，左边的速度u(i,j)为0，系数不使用置为１'
                a_u_add[i, j] = 1

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[i + 1, j] = detar_y * alpha_p/ a_u_add[i + 1, j]
                a_p[i + 1, j] = d[i + 1, j] * detar_y

                '左边界时，左边的速度u(i,j)为0，aw为0'
                d[i - 1, j] = 0
                a_p[i - 1, j] = d[i - 1, j] * detar_y

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]
                a_p[i, j + 1] = d[i, j + 1] * detar_x

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]
                a_p[i, j - 1] = d[i, j - 1] * detar_x

                a_p[i, j] = a_p[i + 1, j] + a_p[i - 1, j] + a_p[i, j + 1] + a_p[i, j - 1]
                b_add = (u[i, j] - u[i + 1, j]) * detar_y + (v[i, j] - v[i, j + 1]) * detar_x

                p_add[i, j] = (a_p[i + 1, j] * p_add[i + 1, j] + a_p[i - 1, j] * p_add[i - 1, j] + a_p[i, j + 1] * p_add[
                                 i, j + 1] + a_p[i, j - 1] * p_add[i, j - 1] + b_add) / a_p[i, j]

                p[i, j] = p[i, j] + alpha_p * p_add[i, j]

            if i == n-1 and j != 0 and j != n - 1 and j != n:
                '右边界时，右边的速度u(i,j)为0，系数不使用置为１'
                a_u_add[i + 1, j] =1

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                '右边界时，右边的速度u(i,j)为0，ae系数为0'
                d[i + 1, j] = 0
                a_p[i + 1, j] = d[i + 1, j] * detar_y

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]
                a_p[i - 1, j] = d[i - 1, j] * detar_y

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]
                a_p[i, j + 1] = d[i, j + 1] * detar_x

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]
                a_p[i, j - 1] = d[i, j - 1] * detar_x

                a_p[i, j] = a_p[i + 1, j] + a_p[i - 1, j] + a_p[i, j + 1] + a_p[i, j - 1]
                b_add = (u[i, j] - u[i + 1, j]) * detar_y + (v[i, j] - v[i, j + 1]) * detar_x

                p_add[i, j] = (a_p[i + 1, j] * p_add[i + 1, j] + a_p[i - 1, j] * p_add[i - 1, j] + a_p[i, j + 1] * p_add[
                                i, j + 1] + a_p[i, j - 1] * p_add[i, j - 1] + b_add) / a_p[i, j]

                p[i, j] = p[i, j] + alpha_p * p_add[i, j]

            if j == 0 and i != 0 and i != n - 1 and i != n:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                '下边界时，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j] = 1

                d[i + 1, j] = detar_y * alpha_p/ a_u_add[i + 1, j]
                a_p[i + 1, j] = d[i + 1, j] * detar_y

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]
                a_p[i - 1, j] = d[i - 1, j] * detar_y

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]
                a_p[i, j + 1] = d[i, j + 1] * detar_x

                '下边界时，下边的速度v(i,j)为0，as为0'
                d[i, j - 1] = 0
                a_p[i, j - 1] = d[i, j - 1] * detar_x

                a_p[i, j] = a_p[i + 1, j] + a_p[i - 1, j] + a_p[i, j + 1] + a_p[i, j - 1]
                b_add = (u[i, j] - u[i + 1, j]) * detar_y + (v[i, j] - v[i, j + 1]) * detar_x

                p_add[i, j] = (a_p[i + 1, j] * p_add[i + 1, j] + a_p[i - 1, j] * p_add[i - 1, j] + a_p[i, j + 1] * p_add[
                                i, j + 1] + a_p[i, j - 1] * p_add[i, j - 1] + b_add) / a_p[i, j]

                p[i, j] = p[i, j] + alpha_p * p_add[i, j]

            if j == n-1 and i != 0 and i != n - 1 and i != n:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                '上边界时，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j + 1] = 1

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[i + 1, j] = detar_y * alpha_p/ a_u_add[i + 1, j]
                a_p[i + 1, j] = d[i + 1, j] * detar_y

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]
                a_p[i - 1, j] = d[i - 1, j] * detar_y

                '上边界时，上边的速度v(i,j)为0，an系数为０'
                d[i, j + 1] = 0
                a_p[i, j + 1] = d[i, j + 1] * detar_x

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]
                a_p[i, j - 1] = d[i, j - 1] * detar_x

                a_p[i, j] = a_p[i + 1, j] + a_p[i - 1, j] + a_p[i, j + 1] + a_p[i, j - 1]
                b_add = (u[i, j] - u[i + 1, j]) * detar_y + (v[i, j] - v[i, j + 1]) * detar_x

                p_add[i, j] = (a_p[i + 1, j] * p_add[i + 1, j] + a_p[i - 1, j] * p_add[i - 1, j] + a_p[i, j + 1] * p_add[
                                  i, j + 1] + a_p[i, j - 1] * p_add[i, j - 1] + b_add) / a_p[i, j]

                p[i, j] = p[i, j] + alpha_p * p_add[i, j]

            '内部节点'
            if i != 0 and i != n - 1 and i != n and j != 0 and j != n - 1 and j != n:
                '补充方程系数'
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re


                d[i+1, j] = detar_y* alpha_u/a_u_add[i+1, j]
                a_p[i+1, j] = d[i+1, j]*detar_y

                d[i-1, j] = detar_y* alpha_u/a_u_add[i, j]
                a_p[i-1, j] = d[i-1, j]*detar_y

                d[i, j+1] = detar_x* alpha_v/a_v_add[i, j+1]
                a_p[i, j+1] = d[i, j+1]*detar_x

                d[i, j-1] = detar_x* alpha_v/a_v_add[i, j]
                a_p[i, j-1] = d[i, j-1]*detar_x

                a_p[i, j] = a_p[i+1, j] + a_p[i-1, j] + a_p[i, j+1] + a_p[i, j-1]
                b_add = (u[i, j] - u[i+1, j])*detar_y + (v[i, j] - v[i, j+1])*detar_x

                p_add[i, j] = (
                              a_p[i + 1, j] * p_add[i + 1, j] + a_p[i - 1, j] * p_add[i - 1, j] + a_p[i, j + 1] * p_add[
                                  i, j + 1] + a_p[i, j - 1] * p_add[i, j - 1] + b_add) / a_p[i, j]

                p[i, j] = p[i, j] + alpha_p * p_add[i, j]



    "------------------------------------------速度修正计算---------------------------------------------"
    for i in range(n + 1):
        for j in range(n + 1):
            '角点'
            if i == 0 and j == 0:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[i, j] = 1

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re
                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j] = 1

                d[i + 1, j] = detar_y * alpha_p/ a_u_add[i + 1, j]

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，aw为0,as=0'
                d[i - 1, j] = 0

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]

                '左下角节点时，左边的速度u(i,j)为0，下边的速度v(i,j)为0，aw为0,as=0'
                d[i, j - 1] = 0

                u[i + 1, j] = u[i + 1, j] + d[i + 1, j]*(p_add[i, j] - p_add[i + 1, j])
                v[i, j + 1] = v[i, j + 1] + d[i, j + 1]*(p_add[i, j] - p_add[i, j + 1])
                u[i - 1, j] = u[i - 1, j] + d[i - 1, j]*(p_add[i - 1, j] - p_add[i, j])
                v[i, j - 1] = v[i, j - 1] + d[i, j - 1]*(p_add[i, j - 1] - p_add[i, j])

            if i == n - 1 and j == 0:
                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[i + 1, j] = 1

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j] = 1

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，ae为0,as为0'
                d[i + 1, j] = 0

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]

                '右下角节点时，右边的速度u(i,j)为0，下边的速度v(i,j)为0，ae为0,as为0'
                d[i, j - 1] = 0

                u[i + 1, j] = u[i + 1, j] + d[i + 1, j]*(p_add[i, j] - p_add[i + 1, j])
                v[i, j + 1] = v[i, j + 1] + d[i, j + 1]*(p_add[i, j] - p_add[i, j + 1])
                u[i - 1, j] = u[i - 1, j] + d[i - 1, j]*(p_add[i - 1, j] - p_add[i, j])
                v[i, j - 1] = v[i, j - 1] + d[i, j - 1]*(p_add[i, j - 1] - p_add[i, j])

            if i == 0 and j == n - 1:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[i, j] = 1

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j + 1] = 1

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[i + 1, j] = detar_y * alpha_p/ a_u_add[i + 1, j]

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，aw为0，an为0'
                d[i - 1, j] = 0

                '左上角节点时，左边的速度u(i,j)为0，上边的速度v(i,j)为0，aw为0，an为0'
                d[i, j + 1] = 0

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]

                u[i + 1, j] = u[i + 1, j] + d[i + 1, j]*(p_add[i, j] - p_add[i + 1, j])
                v[i, j + 1] = v[i, j + 1] + d[i, j + 1]*(p_add[i, j] - p_add[i, j + 1])
                u[i - 1, j] = u[i - 1, j] + d[i - 1, j]*(p_add[i - 1, j] - p_add[i, j])
                v[i, j - 1] = v[i, j - 1] + d[i, j - 1]*(p_add[i, j - 1] - p_add[i, j])

            if i == n - 1 and j == n - 1:
                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_u_add[i + 1, j] = 1

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j + 1] = 1

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，an为0,ae为0'
                d[i + 1, j] = 0

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]

                '右上角节点时，右边的速度u(i,j)为0，上边的速度v(i,j)为0，an为0,ae为0'
                d[i, j + 1] = 0

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]

                u[i + 1, j] = u[i + 1, j] + d[i + 1, j]*(p_add[i, j] - p_add[i + 1, j])
                v[i, j + 1] = v[i, j + 1] + d[i, j + 1]*(p_add[i, j] - p_add[i, j + 1])
                u[i - 1, j] = u[i - 1, j] + d[i - 1, j]*(p_add[i - 1, j] - p_add[i, j])
                v[i, j - 1] = v[i, j - 1] + d[i, j - 1]*(p_add[i, j - 1] - p_add[i, j])

            '边界节点'
            if i == 0 and j!= 0 and j!= n-1 and j!= n:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                '左边界时，左边的速度u(i,j)为0，系数不使用置为１'
                a_u_add[i, j] = 1

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[i + 1, j] = detar_y * alpha_p/ a_u_add[i + 1, j]

                '左边界时，左边的速度u(i,j)为0，aw为0'
                d[i - 1, j] = 0

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]

                u[i + 1, j] = u[i + 1, j] + d[i + 1, j]*(p_add[i, j] - p_add[i + 1, j])
                v[i, j + 1] = v[i, j + 1] + d[i, j + 1]*(p_add[i, j] - p_add[i, j + 1])
                u[i - 1, j] = u[i - 1, j] + d[i - 1, j]*(p_add[i - 1, j] - p_add[i, j])
                v[i, j - 1] = v[i, j - 1] + d[i, j - 1]*(p_add[i, j - 1] - p_add[i, j])

            if i == n-1 and j != 0 and j != n - 1 and j != n:
                '右边界时，右边的速度u(i,j)为0，系数不使用置为１'
                a_u_add[i + 1, j] =1

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                '右边界时，右边的速度u(i,j)为0，ae系数为0'
                d[i + 1, j] = 0

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]

                u[i + 1, j] = u[i + 1, j] + d[i + 1, j]*(p_add[i, j] - p_add[i + 1, j])
                v[i, j + 1] = v[i, j + 1] + d[i, j + 1]*(p_add[i, j] - p_add[i, j + 1])
                u[i - 1, j] = u[i - 1, j] + d[i - 1, j]*(p_add[i - 1, j] - p_add[i, j])
                v[i, j - 1] = v[i, j - 1] + d[i, j - 1]*(p_add[i, j - 1] - p_add[i, j])

            if j == 0 and i != 0 and i != n - 1 and i != n:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                '下边界时，下边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j] = 1

                d[i + 1, j] = detar_y * alpha_p/ a_u_add[i + 1, j]

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]

                d[i, j + 1] = detar_x * alpha_p/ a_v_add[i, j + 1]

                '下边界时，下边的速度v(i,j)为0，as为0'
                d[i, j - 1] = 0

                u[i + 1, j] = u[i + 1, j] + d[i + 1, j]*(p_add[i, j] - p_add[i + 1, j])
                v[i, j + 1] = v[i, j + 1] + d[i, j + 1]*(p_add[i, j] - p_add[i, j + 1])
                u[i - 1, j] = u[i - 1, j] + d[i - 1, j]*(p_add[i - 1, j] - p_add[i, j])
                v[i, j - 1] = v[i, j - 1] + d[i, j - 1]*(p_add[i, j - 1] - p_add[i, j])

            if j == n-1 and i != 0 and i != n - 1 and i != n:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                '上边界时，上边的速度v(i,j)为0，系数不使用置为１'
                a_v_add[i, j + 1] = 1

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[i + 1, j] = detar_y * alpha_p/ a_u_add[i + 1, j]

                d[i - 1, j] = detar_y * alpha_p/ a_u_add[i, j]

                '上边界时，上边的速度v(i,j)为0，an系数为０'
                d[i, j + 1] = 0

                d[i, j - 1] = detar_x * alpha_p/ a_v_add[i, j]

                u[i + 1, j] = u[i + 1, j] + d[i + 1, j]*(p_add[i, j] - p_add[i + 1, j])
                v[i, j + 1] = v[i, j + 1] + d[i, j + 1]*(p_add[i, j] - p_add[i, j + 1])
                u[i - 1, j] = u[i - 1, j] + d[i - 1, j]*(p_add[i - 1, j] - p_add[i, j])
                v[i, j - 1] = v[i, j - 1] + d[i, j - 1]*(p_add[i, j - 1] - p_add[i, j])


            '内部节点'
            if i != 0 and i != n - 1 and i != n and j != 0 and j != n - 1 and j != n:
                a_u_add[i + 1, j] = (vv[i + 1, j] + vv[i, j] + vv[i, j + 1] + vv[i + 1, j + 1]) / 4 * detar_x \
                                    + uu[i + 1, j] * detar_y \
                                    + detar_x * detar_y / detar_t \
                                    + detar_y / detar_x * 2 / Re \
                                    + detar_x / detar_y * 2 / Re

                a_u_add[i, j] = (vv[i, j] + vv[i - 1, j] + vv[i - 1, j + 1] + vv[i, j + 1]) / 4 * detar_x \
                                + uu[i, j] * detar_y \
                                + detar_x * detar_y / detar_t \
                                + detar_y / detar_x * 2 / Re \
                                + detar_x / detar_y * 2 / Re

                a_v_add[i, j + 1] = (uu[i + 1, j + 1] + uu[i, j + 1] + uu[i, j] + uu[i + 1, j]) / 4 * detar_y * detar_t \
                                    + vv[i, j + 1] * detar_x * detar_t \
                                    + detar_x * detar_y \
                                    + detar_y * detar_t / detar_x * 2 / Re \
                                    + detar_x * detar_t / detar_y * 2 / Re

                a_v_add[i, j] = (uu[i + 1, j] + uu[i, j] + uu[i, j - 1] + uu[i + 1, j - 1]) / 4 * detar_y * detar_t \
                                + vv[i, j] * detar_x * detar_t \
                                + detar_x * detar_y \
                                + detar_y * detar_t / detar_x * 2 / Re \
                                + detar_x * detar_t / detar_y * 2 / Re

                d[i + 1, j] = detar_y* alpha_p/a_u_add[i + 1, j]

                d[i, j + 1]= detar_x* alpha_p/a_v_add[i, j + 1]

                d[i - 1, j]= detar_y* alpha_p/a_u_add[i, j]

                d[i, j - 1]= detar_x* alpha_p/a_v_add[i, j]

                u[i + 1, j] = u[i + 1, j] + d[i + 1, j]*(p_add[i, j] - p_add[i + 1, j])
                v[i, j + 1] = v[i, j + 1] + d[i, j + 1]*(p_add[i, j] - p_add[i, j + 1])
                u[i - 1, j] = u[i - 1, j] + d[i - 1, j]*(p_add[i - 1, j] - p_add[i, j])
                v[i, j - 1] = v[i, j - 1] + d[i, j - 1]*(p_add[i, j - 1] - p_add[i, j])
    '将本轮迭代的u,v保存'
    vv = v
    uu = u
# print(p)
print(u)
print(v)
# print(a_u)
# print(a_v)
# print(a_u_add)
# print(a_v_add)
# print(a_p)
# print(p_add)