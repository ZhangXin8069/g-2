import numpy as np
import matplotlib.pyplot as plt

# 设置问题区域和网格参数
x_start, x_end = 0, 2
y_start, y_end = x_end + 0.01, 30 + 0.01

x_resolution, y_resolution = 1000, 100

hx = (x_end - x_start) / (x_resolution - 1)
hy = (y_end - y_start) / (y_resolution - 1)

x_grid = np.linspace(x_start, x_end, x_resolution)
y_grid = np.linspace(y_start, y_end, y_resolution)

# 定义函数f


def f(x, M, G, F):
    return F**2 * ((M**4 + M**2 * G**2) / (24 * np.pi * ((x - M**2)**2 + M**2 * G**2)))


# 参数设置
params = [
    (0.22, 0.78, 0.15),
    (0.19, 1.46, 0.4),
    (0.14, 1.72, 0.25),
    (0.14, 1.90, 0.1)
]

# 计算F
F = np.sum([f(x_grid, M, G, FV) for FV, M, G in params], axis=0)

# 计算G1和G2
G1 = np.trapz(np.sin(2 * np.pi * x_grid) /
              (y_grid[:, None] - x_grid), x_grid, axis=1)
G2 = np.trapz((x_grid - x_grid) / (y_grid[:, None] - x_grid), x_grid, axis=1)

# 定义和计算误差参数和精确解
B1, B2 = 1, 2
error_scale = 1.01

b1, b2 = error_scale * B1, error_scale * B2

f_exact = B1 * np.sin(2 * np.pi * x_grid) + B2 * (x_grid - x_grid)
f_delta = b1 * np.sin(2 * np.pi * x_grid) + b2 * (x_grid - x_grid)

G1_exact, G2_exact = B1 * G1, B2 * G2
G1_delta, G2_delta = b1 * G1, b2 * G2

# 计算V_exact和G3
V_exact = (f_delta[0] * (x_grid - x_end) / (x_start - x_end)) + \
    (f_delta[-1] * (x_grid - x_start) / (x_end - x_start))
G3 = np.trapz(V_exact / (y_grid[:, None] - x_grid), x_grid, axis=1)

# 计算U_exact和G_exact
U_exact = f_exact - V_exact
G_exact = G1_exact + G2_exact - G3

# 构建矩阵A和H1
Aphi = np.zeros((y_resolution, x_resolution))

for j in range(1, x_resolution - 1):
    Aphi[:, j] = (1 / hx) * ((y_grid - x_grid[j-1]) * np.log(y_grid - x_grid[j-1]) +
                             (x_grid[j-1] - x_grid[j]) + (y_grid - x_grid[j+1]) * np.log(y_grid - x_grid[j+1]) -
                             2 * (y_grid - x_grid[j]) * np.log(y_grid - x_grid[j]) - (x_grid[j] - x_grid[j+1]))

Aphi[:, 0] = np.log(y_grid - x_grid[0]) + (1 / hx) * ((y_grid - x_grid[1]) * np.log(y_grid - x_grid[1]) -
                                                      (y_grid - x_grid[0]) * np.log(y_grid - x_grid[0]) -
                                                      (x_grid[0] - x_grid[1]))
Aphi[:, -1] = -np.log(y_grid - x_grid[-1]) + (1 / hx) * (-(y_grid - x_grid[-1]) * np.log(y_grid - x_grid[-1]) +
                                                         (y_grid - x_grid[-2]) * np.log(y_grid - x_grid[-2]) +
                                                         (x_grid[-2] - x_grid[-1]))

A = np.zeros((x_resolution - 2, x_resolution - 2))
ba = np.zeros((x_resolution - 2, 1))

M = np.zeros((x_resolution - 2, x_resolution - 2))
M1 = np.zeros((x_resolution - 2, x_resolution - 2))

for j in range(1, x_resolution - 1):
    A[:, j-1] = np.trapz(Aphi[:, j] * Aphi[:, j-1], y_grid)

for i in range(1, x_resolution - 1):
    for j in range(1, x_resolution - 1):
        M[i-1, j-1] = np.trapz(Aphi[:, i] * Aphi[:, j], y_grid)
        M1[i-1, j-1] = 2 * np.trapz(Aphi[:, i] * Aphi[:, j], y_grid)

# 构建矩阵H1
H1 = (hx / 12) * (M + M1)

# 删除积分矩阵中的第一列和最后一列对应的元素，因为它们都为0
A = np.delete(A, [0, -1], axis=1)
H1 = np.delete(H1, [0, -1], axis=0)
H1 = np.delete(H1, [0, -1], axis=1)

# 迭代选择正则化参数
num_iterations = 20
alpha_values = np.zeros(num_iterations)
L_curve_values = np.zeros(num_iterations)
LH_curve_values = np.zeros(num_iterations)
GCV_values = np.zeros(num_iterations)
Rho_values = np.zeros(num_iterations)

for k in range(num_iterations):
    alpha_values[k] = 10**(-k)  # 选择不同的正则化参数 alpha

    ua = np.linalg.solve(A + alpha_values[k] * H1, ba)
    Ua = np.concatenate(([0], ua, [0]))
    fa = Ua + V_exact

    # 计算误差和指标
    res = np.sqrt(np.trapz((Aphi @ fa - G1_delta - G2_delta)**2, y_grid))
    err = np.sqrt(np.trapz((fa - f_exact)**2, x_grid))

    # 计算 L 曲线和 LH 曲线的值
    L_curve_values[k] = np.sqrt(np.trapz(fa**2, x_grid)) * res
    LH_curve_values[k] = np.sqrt(
        np.trapz(fa**2 + np.gradient(fa)**2, x_grid)) * res

    # 使用 GCV 方法选择最佳 alpha
    VV = np.eye(y_resolution) - Aphi @ np.linalg.solve(Aphi.T @
                                                       Aphi + alpha_values[k] * np.eye(x_resolution), Aphi.T)
    V1 = np.linalg.norm(VV @ (G1_delta + G2_delta))**2
    V2 = np.trace(VV)**2
    GCV_values[k] = V1 / V2

    # 计算 Rho
    Rho_values[k] = np.linalg.norm(alpha_values[k] * np.linalg.solve(
        Aphi.T @ Aphi + alpha_values[k] * np.eye(x_resolution), fa))

# 选择最佳的 alpha 值和方法
min_L_curve_index = np.argmin(L_curve_values)
min_LH_curve_index = np.argmin(LH_curve_values)
min_GCV_index = np.argmin(GCV_values)
min_Rho_index = np.argmin(Rho_values)

# 输出结果
result_methods = ['L-curve', 'LH-curve', 'GCV', '拟最优']
result_alpha_values = [alpha_values[min_L_curve_index], alpha_values[min_LH_curve_index],
                       alpha_values[min_GCV_index], alpha_values[min_Rho_index]]

results = {
    '选择方法': result_methods,
    '最佳 alpha 值': result_alpha_values
}

print(results)  # 显示结果
plt.show()  # 显示图像
