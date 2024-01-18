import numpy as np
import intvalpy as ip
import matplotlib.pyplot as plt
from intvalpy.utils import asinterval, intersection, dist, infinity, isnan
from intvalpy.RealInterval import ARITHMETICS

eps1 = 0.05
eps2 = 0.1

def f1(x, y):
    return x**2 + y**2 - 1 
def f2(x, y):
    return x-y**2 + 0.5 - ip.Interval(0, eps2, midRadQ=True)
def df1_x(x):
    return 2*x
def df1_y(y):
    return 2*y
def df2_x(x):
    return 1
def df2_y(y):
    return -2*y

def J(X):
    midL = np.zeros((2, 2))
    radL = np.zeros((2, 2))
    
    radL[0][0] = df1_x(X[0]).rad
    midL[0][0] = df1_x(X[0]).mid

    radL[0][1] = df1_y(X[1]).rad
    midL[0][1] = df1_y(X[1]).mid

    midL[1][0] = 1
    radL[1][0] = 0

    radL[1][1] = df2_y(X[1]).rad
    midL[1][1] = df2_y(X[1]).mid
    return ip.Interval(midL, radL, midRadQ=True)

def F(X):
    ere = f1(X[0],X[1])
    return [f1(X[0],X[1]), f2(X[0],X[1])]

def spectral_radius(A):
    eigenvalues = np.linalg.eigvals(A)
    spectral_radius = np.max(np.abs(eigenvalues))
    return spectral_radius

def Krawczyk(func, J, x0, maxiter=2000, tol=1e-5):
    spect_rad =[]
    def K(X, c):
        L = J(X)
        spect_rad.append(spectral_radius(L).mid)
        LAMBDA = np.linalg.inv(L.to_float().mid)
        B = np.eye(2) - LAMBDA @ L
        w, _ = np.linalg.eigh(B.to_float().mag)
        return c - LAMBDA @ func(c) + B @ (X - c)

    result = x0
    pre_result = result.copy
    c = asinterval(result.mid)

    error = infinity
    nit = 0
    X_k = []
    X_k.append(result)
    Kr = []
    i = 0
    while nit <= maxiter and error > tol:
        i += 1
        krav = K(result, c)
        Kr.append(krav)
        result = intersection(result, krav)

        if isnan(result).any():
            X_k.append(result)
            return result
        X_k.append(result)
        c = asinterval(result.mid)
        error = dist(result, pre_result)
        pre_result = result.copy
        nit += 1

    return result, X_k, Kr, i,spect_rad

def print_int(X):
    for i in range(len(X)):
        a_1 = float(X_k[i][0].a)
        a_2 = float(X_k[i][1].a)
        b_1 = float(X_k[i][0].b)
        b_2 = float(X_k[i][1].b)
        print(i+1, " &$ [", "{:.6}".format(a_1),",", "{:.6}".format(b_1), "]$ & $[", "{:.6}".format(a_2), ",", "{:.6}".format(b_2), "]$ &" ,"{:.6}".format(b_1-a_1), "&" "{:.6}".format(b_2 - a_2), "&\\\\" )

if __name__ == "__main__":

    midX = [0.5, 0.5]
    radX = [0.5, 0.5]

    X = ip.Interval(midX, radX, midRadQ=True)
    ttt = ip.Interval(0, 0.2, midRadQ=True)

    z, X_k, Kr, i, spect_rad = Krawczyk(F, J, X)
    
    for i in range(len(X_k)):
        one = abs(X_k[i][0].b - X_k[i][0].a)
        two = abs(X_k[i][1].b - X_k[i][1].a)
        iveRect = plt.Rectangle((X_k[i][0].a, X_k[i][1].a), one , two, edgecolor='black', facecolor='none', label='Брус ive', linewidth=1.5)
        plt.gca().add_patch(iveRect)



    print("Кол-во итераций", i)
    plt.plot(midX[0], midX[1], 'r*', ms=4, label='start')
    print_int(X_k)
    x = np.arange(-2, 2, 0.01)
    y1 = np.sqrt(1 - pow(x, 2) )
    plt.plot(x, y1, '--', linewidth=0.7)
    y2 = np.sqrt(1 - pow(x, 2) )
    plt.plot(x, y2, '--', linewidth=0.7)
    y3 = np.sqrt(x + 0.5 - eps2 )
    plt.plot(x, y3, linewidth=0.7)
    y4 = np.sqrt(x + 0.5 + eps2 )
    plt.plot(x, y4, linewidth=0.7)
    plt.grid()
    plt.xlim(0, 2)
    plt.ylim(0, 2)
    plt.show()
    print(spect_rad)