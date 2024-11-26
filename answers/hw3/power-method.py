import numpy as np

np.random.seed(20241025) # Set random seed.
A = np.random.randn(2000, 1000)


def Power_method(A, k, Iterations=100):
    m, n = A.shape
    # print("A.shape: ", A.shape)
    x = np.random.randn(n)
    # print("x.shape: ", x.shape)
    Q = np.zeros((n, k))
    A_T = A.T
    B = A_T @ A
    for i in range(k):
        Q[:, i] = x / np.linalg.norm(x)
        x = B @ x
    for i in range(Iterations):
        Q = B @ Q
        # 对Q正交化
        Q, _ = np.linalg.qr(Q)
    singular_values = np.zeros(k)
    for i in range(k):
        singular_values[i] = np.linalg.norm(A @ Q[:, i])
    # 逆向排序
    singular_values = np.sort(singular_values)[::-1]
    return singular_values


k = 10
i = 10000
singular_values = Power_method(A, k, i)
print ("迭代次数：", i)
print(f"幂方法前{k}个奇异值：") # 精确到两位小数
print(np.round(singular_values, 2))
singular_values = np.linalg.svd(A, compute_uv=False)
print(f"svd调包前{k}个奇异值：")
print(np.round(singular_values[:k],2))
