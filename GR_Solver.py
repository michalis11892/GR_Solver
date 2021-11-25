import sympy as sym
import numpy as np
import copy 

#Defining Coordinates
x = [sym.symbols('t'), sym.symbols('r'), sym.symbols('theta'), sym.symbols('phi')]

#Defining Constants
G, M = sym.symbols('G, M')

#Defining the Metric Tensor (g_i_j)
g = [[-1*(1-2*G*M/x[1]), 0, 0, 0],
     [0, 1/(1-2*G*M/x[1]), 0, 0],
     [0, 0, x[1]**2, 0],
     [0, 0, 0, x[1]**2 *sym.sin(x[2])**2]]

print('Metric (g_i_j): \n', np.array(g), '\n')

#Calculating the Inverse Metric Tensor (g^i^j)
g_inv = sym.Matrix(g).inv()

print('Inverse Metric (g^i^j): \n', np.array(g_inv), '\n')

#Defining the Christofel Pseudo-Tensor (gamma^i_j_k)
gamma = []
for i in range(4):
    m_temp = []
    for j in range(4):
        temp = []
        for k in range(4):
            temp.append(0)
        m_temp.append(temp)
    gamma.append(m_temp)

#Calculating the Christofel Pseudo-Tensor (gamma^i_j_k)
for i in range(4):
    for j in range(4):
        for k in range(4):
            for s in range(4):
                gamma[i][j][k] += 0.5*g_inv[i, s]*(sym.diff(g[k][s], x[j])+sym.diff(g[j][s], x[k])-sym.diff(g[j][k], x[s]))
            gamma[i][j][k] = sym.simplify(gamma[i][j][k]) #Expression Simplification

print('Christofel Pseudo-Tensor (gamma^i_j_k): \n', np.array(gamma), '\n')

#Defining the Riemann Curvature Tensor (R^i_j_k_w)
R = []
for i in range(4):
    b_temp = []
    for j in range(4):
        m_temp = []
        for k in range(4):
            temp = []
            for w in range(4):
                temp.append(0)
            m_temp.append(temp)
        b_temp.append(m_temp)
    R.append(b_temp)

#Defining the Riemann Curvature Tensor (R_i_j_k_w)
R_ = copy.deepcopy(R)

#Defining the Riemann Curvature Tensor (R^i^j^k^w)
Ru = copy.deepcopy(R) 

#Calculating the Riemann Curvature Tensor (R^i_j_k_w)
for i in range(4):
    for j in range(4):
        for k in range(4):
            for w in range(4):
                R[i][j][k][w] = sym.diff(gamma[i][j][k], x[w]) - sym.diff(gamma[i][j][w], x[k])
                for s in range(4):
                    R[i][j][k][w] += gamma[s][j][k]*gamma[i][w][s] - gamma[s][j][w]*gamma[i][k][s]
                R[i][j][k][w] = sym.simplify(R[i][j][k][w]) #Expression Simplification

print('Riemann Curvature Tensor (R^i_j_k_w): \n', np.array(R), '\n')

#Calculating the Riemann Curvature Tensor (R_i_j_k_w)
for i in range(4):
    for j in range(4):
        for k in range(4):
            for w in range(4):
                for s in range(4):
                    R_[i][j][k][w] += g[i][s]*R[s][j][k][w]
                R_[i][j][k][w] = sym.simplify(R_[i][j][k][w]) #Expression Simplification

print('Riemann Curvature Tensor (R_i_j_k_w): \n', np.array(R_), '\n')

#Calculating the Riemann Curvature Tensor (R^i^j^k^w)
for i in range(4):
    for j in range(4):
        for k in range(4):
            for w in range(4):
                for s1 in range(4):
                    for s2 in range(4):
                        for s3 in range(4):
                            Ru[i][j][k][w] += g_inv[j, s1]*g_inv[k, s2]*g_inv[w, s3]*R[i][s1][s2][s3]
                Ru[i][j][k][w] = sym.simplify(Ru[i][j][k][w]) #Expression Simplification

print('Riemann Curvature Tensor (R^i^j^k^w): \n', np.array(Ru), '\n')

#Defining the Ricci Tensor (R_i_j)
r = []
for i in range(4):
    temp = []
    for j in range(4):
        temp.append(0)
    r.append(temp)

#Calculating the Ricci Tensor (R_i_j)
for i in range(4):
    for j in range(4):
        for s1 in range(4):
            for s2 in range(4):
                r[i][j] += g_inv[s1, s2]*R_[s1][i][s2][j]
            r[i][j] = sym.simplify(r[i][j]) #Expression Simplification


print('Ricci Tensor (R_i_j): \n', np.array(r), '\n')

#Defining & Calculating the Scalar Curvature (R)
K = 0
for s1 in range(4):
    for s2 in range(4):
        K += g_inv[s1, s2]*r[s1][s2]
K = sym.simplify(K) #Expression Simplification

print('Scalar Curvature (R): \n', K, '\n')

#Defining & Calculating the Kretschmann Scalar (R_i_j_k_w*R^i^j^k^w):
RR = 0
for s1 in range(4):
    for s2 in range(4):
        for s3 in range(4):
            for s4 in range(4):
                RR += R_[s1][s2][s3][s4]*Ru[s1][s2][s3][s4]
RR = sym.simplify(RR) #Expression Simplification

print('Kretschmann Scalar (R_i_j_k_w*R^i^j^k^w): \n', RR, '\n')
