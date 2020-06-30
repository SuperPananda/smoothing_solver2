import numpy as np
import seaborn as sns 
import pandas as pd
from scipy import sparse
from math import sqrt, fabs, log
import matplotlib.pyplot as plt
fig = plt.subplots()

""" def createglobal(matrix,i1,i2):
    for i in range(4):
        for j in range(4):
            glmatrix[i1] """
#r, z, temp, L, U, LUdi, tt = np.array([])
maxiter = 1000
eps = 1.0E-18


def begin():
    global r, z, temp, tt
    r = np.zeros(len(Setpoints))
    z = np.zeros(len(Setpoints))
    temp = np.zeros(len(Setpoints))
    tt = np.zeros(len(Setpoints))


def mult(a, res):
    for i in range(len(Setpoints)):
        i0 = ia[i]
        i1 = ia[i+1]
        res[i] = di[i]*a[i]
        for k in range(i0, i1):
            j = ja[k]
            res[i] += gg[k]*a[j]
            res[j] += gg[k]*a[i]
    return res


def scalar_mult(a, b):
    s = 0.0
    for i in range(len(Setpoints)):
        s += a[i]*b[i]
    return s


def norma(vect):
    return sqrt(scalar_mult(vect, vect))


def calc_nev(res):
    res = mult(x, res)
    for i in range(len(Setpoints)):
        res[i] = globalf[i]-res[i]
    return res


def create_LU():
    global L, U, LUdi
    begin()
    size = ia[len(Setpoints)]
    L = np.zeros(size)
    U = np.zeros(size)
    LUdi = np.zeros(len(Setpoints))
    for i in range(len(Setpoints)):
        sd = 0
        i0 = ia[i]
        i1 = ia[i+1]
        for k in range(i0, i1):
            su = 0
            sl = 0
            j = ja[k]
            kj = ia[j]
            ki = ia[i]
            j1 = ia[j + 1]
            while ki < k and kj < j1:
                if ja[kj] == ja[ki]:
                    sl += L[ki] * U[kj]
                    su += L[kj] * U[ki]
                    ki = ki+1
                    kj = kj+1
                else:
                    if ja[ki] < ja[kj]:
                        ki = ki+1
                    else:
                        kj = kj+1
            L[k] = gg[k] - sl
            U[k] = (gg[k] - su) / LUdi[j]
            sd += L[k] * U[k]
        LUdi[i] = di[i] - sd


def mult_U1(y):
    for i in reversed(range(len(Setpoints))):
        i0 = ia[i]
        i1 = ia[i + 1]
        xi = y[i]
        for k in range(i0, i1):
            y[ja[k]] -= U[k] * xi
    return y


def mult_L1(y):
    for i in range(len(Setpoints)):
        i0 = ia[i]
        i1 = ia[i + 1]
        for k in range(i0, i1):
            y[i] -= y[ja[k]] * L[k]
        y[i] /= LUdi[i]
    return y


def los_with_LU():
    global x, r, z, tt, temp
    p = np.zeros(len(Setpoints))
    create_LU()
    r = calc_nev(r)
    r = mult_L1(r)
    for i in range(len(Setpoints)):
        z[i] = r[i]
    z = mult_U1(z)
    p = mult(z, p)
    p = mult_L1(p)
    norm_f = norma(r)
    nev = norma(r) / norm_f
    k = 1
    while nev > eps and k < maxiter:
        p_na_p = scalar_mult(p, p)
        alfa = scalar_mult(p, r) / p_na_p
        for i in range(len(Setpoints)):
            x[i] += alfa * z[i]
            r[i] -= alfa * p[i]
        for i in range(len(Setpoints)):
            tt[i] = r[i]
        tt = mult_U1(tt)
        temp = mult(tt, temp)
        temp = mult_L1(temp)
        betta = -scalar_mult(p, temp) / p_na_p
        for i in range(len(Setpoints)):
            p[i] = temp[i] + betta * p[i]
            z[i] = tt[i] + betta * z[i]
        nev = norma(r)/norm_f


def sborka(femel, a, b):
    global di, globalf, gg
    for i in range(16):
        ki = femel[i]
        di[ki] += a[i][i]
        globalf[ki] += b[i]
        for j in range(i):
            kj = femel[j]
            if ki > kj:
                k = ia[ki]
            else:
                k = ia[kj]
                kj = ki
            while ja[k] != kj:
                k = k+1
            gg[k] += a[i][j]


def add_str(femel):
    global ig
    for i in range(16):
        for j in range(16):
            ki = femel[i]
            kj = femel[j]
            # t = kj in ig[ki]
            if ki > kj and ((kj in ig[ki]) == False):
                ig[ki].append(kj)


def create_ia_ja(NumFiniteElem, NumCoord):
    global ia, gg, ja
    for i in range(NumFiniteElem):
        add_str(finite[i])

    ia = np.zeros(NumCoord+1)
    ia = np.int_(ia)
    ia[0] = 0

    for i in range(1, NumCoord+1):
        ia[i] = ia[i-1] + len(ig[i-1])

    ja = np.zeros(ia[NumCoord])
    ja = np.int_(ja)

    i = 0
    for j in range(NumCoord):
        # t = not ig[j]
        while ((not ig[j]) == False):
            ja[i] = ig[j].pop(0)
            i = i+1
    gg = np.zeros(len(ja))


""" def d2phi(i,x):
    if i==0:
        fi=(-1)*(6*(h-2*x+2*x1))/(h*h*h)
    elif i==1:
        fi = -4
    elif i==2:
        fi = 3*xsi*xsi+2*xsi*xsi*xsi
    elif i==3:
        fi = 2*((-1)*xsi*xsi+xsi*xsi*xsi)
    return fi """


def phi(i, xsi):
    """ xsi = (x-2)/2 """
    if i == 0:
        fi = -9/2*(xsi-1/3)*(xsi-2/3)*(xsi-1)
    elif i == 1:
        fi = 27/2*xsi*(xsi-2/3)*(xsi-1)
    elif i == 2:
        fi = -27/2*xsi*(xsi-1/3)*(xsi-1)
    elif i == 3:
        fi = 9/2*(xsi-1/3)*(xsi-2/3)*xsi
    return fi


def function_y(i):
    return (((i+1)-1) % 4)+1-1


def function_v(i):
    return (((i+1)-1)//4)+1-1


def b_func(i, x, x1, x2):
    if i == 0:
        fi = (x2 - x)/(x2-x1)
    elif i == 1:
        fi = (x-x1)/(x2-x1)
    return fi

def b_dfunc(i,x,x1,x2):
    h = x2-x1
    if i==0:
        fi = -1/h
    elif i==1:
        fi = 1/h
    return fi

def d_pfi_dx(i,x,x1,x2,y,y1,y2):
    if i == 0:
        fi = b_dfunc(0, x, x1, x2)*b_func(0, y, y1, y2)
    elif i == 1:
        fi = b_dfunc(1, x, x1, x2)*b_func(0, y, y1, y2)
    elif i == 2:
        fi = b_dfunc(0, x, x1, x2)*b_func(1, y, y1, y2)
    elif i == 3:
        fi = b_dfunc(1, x, x1, x2)*b_func(1, y, y1, y2)
    return fi

def d_pfi_dy(i,x,x1,x2,y,y1,y2):
    if i == 0:
        fi = b_func(0, x, x1, x2)*b_dfunc(0, y, y1, y2)
    elif i == 1:
        fi = b_func(1, x, x1, x2)*b_dfunc(0, y, y1, y2)
    elif i == 2:
        fi = b_func(0, x, x1, x2)*b_dfunc(1, y, y1, y2)
    elif i == 3:
        fi = b_func(1, x, x1, x2)*b_dfunc(1, y, y1, y2)
    return fi

def d_pfi_lag_dx(i,x,x1,x2,y,y1,y2):
    if i == 0:
        fi = db_func_lag(0,x,x1,x2)*b_func_lag(0,y,y1,y2)
    elif i == 1:
        fi = db_func_lag(1,x,x1,x2)*b_func_lag(0,y,y1,y2)
    elif i == 2:
        fi = db_func_lag(2,x,x1,x2)*b_func_lag(0,y,y1,y2)
    elif i == 3:
        fi = db_func_lag(3,x,x1,x2)*b_func_lag(0,y,y1,y2)
    elif i == 4:
        fi = db_func_lag(0,x,x1,x2)*b_func_lag(1,y,y1,y2)
    elif i == 5:
        fi = db_func_lag(1,x,x1,x2)*b_func_lag(1,y,y1,y2)
    elif i == 6:
        fi = db_func_lag(2,x,x1,x2)*b_func_lag(1,y,y1,y2)
    elif i == 7:
        fi = db_func_lag(3,x,x1,x2)*b_func_lag(1,y,y1,y2)
    elif i == 8:
        fi = db_func_lag(0,x,x1,x2)*b_func_lag(2,y,y1,y2)
    elif i == 9:
        fi = db_func_lag(1,x,x1,x2)*b_func_lag(2,y,y1,y2)
    elif i == 10:
        fi = db_func_lag(2,x,x1,x2)*b_func_lag(2,y,y1,y2)
    elif i == 11:
        fi = db_func_lag(3,x,x1,x2)*b_func_lag(2,y,y1,y2)
    elif i == 12:
        fi = db_func_lag(0,x,x1,x2)*b_func_lag(3,y,y1,y2)
    elif i == 13:
        fi = db_func_lag(1,x,x1,x2)*b_func_lag(3,y,y1,y2)
    elif i == 14:
        fi = db_func_lag(2,x,x1,x2)*b_func_lag(3,y,y1,y2)
    elif i == 15:
        fi = db_func_lag(3,x,x1,x2)*b_func_lag(3,y,y1,y2)
    return fi

def d_pfi_lag_dy(i,x,x1,x2,y,y1,y2):
    if i == 0:
        fi = b_func_lag(0,x,x1,x2)*db_func_lag(0,y,y1,y2)
    elif i == 1:
        fi = b_func_lag(1,x,x1,x2)*db_func_lag(0,y,y1,y2)
    elif i == 2:
        fi = b_func_lag(2,x,x1,x2)*db_func_lag(0,y,y1,y2)
    elif i == 3:
        fi = b_func_lag(3,x,x1,x2)*db_func_lag(0,y,y1,y2)
    elif i == 4:
        fi = b_func_lag(0,x,x1,x2)*db_func_lag(1,y,y1,y2)
    elif i == 5:
        fi = b_func_lag(1,x,x1,x2)*db_func_lag(1,y,y1,y2)
    elif i == 6:
        fi = b_func_lag(2,x,x1,x2)*db_func_lag(1,y,y1,y2)
    elif i == 7:
        fi = b_func_lag(3,x,x1,x2)*db_func_lag(1,y,y1,y2)
    elif i == 8:
        fi = b_func_lag(0,x,x1,x2)*db_func_lag(2,y,y1,y2)
    elif i == 9:
        fi = b_func_lag(1,x,x1,x2)*db_func_lag(2,y,y1,y2)
    elif i == 10:
        fi = b_func_lag(2,x,x1,x2)*db_func_lag(2,y,y1,y2)
    elif i == 11:
        fi = b_func_lag(3,x,x1,x2)*db_func_lag(2,y,y1,y2)
    elif i == 12:
        fi = b_func_lag(0,x,x1,x2)*db_func_lag(3,y,y1,y2)
    elif i == 13:
        fi = b_func_lag(1,x,x1,x2)*db_func_lag(3,y,y1,y2)
    elif i == 14:
        fi = b_func_lag(2,x,x1,x2)*db_func_lag(3,y,y1,y2)
    elif i == 15:
        fi = b_func_lag(3,x,x1,x2)*db_func_lag(3,y,y1,y2)
    return fi

def b_func_lag(i,x,x1,x2):
    xsi = (x-x1)/(x2-x1)
    if i == 0:
        fi = -9/2 * (xsi-1/3)*(xsi-2/3)*(xsi-1)
    elif i == 1:
        fi = 27/2 * xsi*(xsi-2/3)*(xsi-1)
    elif i == 2:
        fi = -27/2 * xsi*(xsi-1/3)*(xsi-1)
    elif i == 3:
        fi = 9/2 * (xsi-1/3)*(xsi-2/3)*xsi
    return fi

def db_func_lag(i,x,x1,x2):
    h = x2-x1
    xsi = (x-x1)/h
    if i==0:
        fi = 1/2*(-27*xsi*xsi+36*xsi-11)*1/h
    elif i==1:
        fi = (81*xsi*xsi/2-45*xsi+9)*1/h
    elif i==2:
        fi = -9/2*(9*xsi*xsi-8*xsi+1)*1/h
    elif i == 3:
        fi = (27*xsi*xsi/2-9*xsi+1)*1/h
    return fi

def pfi_lag(i,x,x1,x2,y,y1,y2):
    if i == 0:
        fi = b_func_lag(0,x,x1,x2)*b_func_lag(0,y,y1,y2)
    elif i == 1:
        fi = b_func_lag(1,x,x1,x2)*b_func_lag(0,y,y1,y2)
    elif i == 2:
        fi = b_func_lag(2,x,x1,x2)*b_func_lag(0,y,y1,y2)
    elif i == 3:
        fi = b_func_lag(3,x,x1,x2)*b_func_lag(0,y,y1,y2)
    elif i == 4:
        fi = b_func_lag(0,x,x1,x2)*b_func_lag(1,y,y1,y2)
    elif i == 5:
        fi = b_func_lag(1,x,x1,x2)*b_func_lag(1,y,y1,y2)
    elif i == 6:
        fi = b_func_lag(2,x,x1,x2)*b_func_lag(1,y,y1,y2)
    elif i == 7:
        fi = b_func_lag(3,x,x1,x2)*b_func_lag(1,y,y1,y2)
    elif i == 8:
        fi = b_func_lag(0,x,x1,x2)*b_func_lag(2,y,y1,y2)
    elif i == 9:
        fi = b_func_lag(1,x,x1,x2)*b_func_lag(2,y,y1,y2)
    elif i == 10:
        fi = b_func_lag(2,x,x1,x2)*b_func_lag(2,y,y1,y2)
    elif i == 11:
        fi = b_func_lag(3,x,x1,x2)*b_func_lag(2,y,y1,y2)
    elif i == 12:
        fi = b_func_lag(0,x,x1,x2)*b_func_lag(3,y,y1,y2)
    elif i == 13:
        fi = b_func_lag(1,x,x1,x2)*b_func_lag(3,y,y1,y2)
    elif i == 14:
        fi = b_func_lag(2,x,x1,x2)*b_func_lag(3,y,y1,y2)
    elif i == 15:
        fi = b_func_lag(3,x,x1,x2)*b_func_lag(3,y,y1,y2)
    return fi

def pfi(i, x, x1, x2, y, y1, y2):
    # phi(function_y(i),xsi)*phi(function_y(i),nsi)
    if i == 0:
        fi = b_func(0, x, x1, x2)*b_func(0, y, y1, y2)
    elif i == 1:
        fi = b_func(1, x, x1, x2)*b_func(0, y, y1, y2)
    elif i == 2:
        fi = b_func(0, x, x1, x2)*b_func(1, y, y1, y2)
    elif i == 3:
        fi = b_func(1, x, x1, x2)*b_func(1, y, y1, y2)
    return fi


def createmesh():
    global finite, Setpoints
    kx = 9
    ky = 9
    coord = np.loadtxt('in.txt')
    x0 = coord[0][0]
    y0 = coord[0][1]

    hx = (coord[0][2] - coord[0][0]) / (4*kx-kx)
    hy = (coord[1][1] - coord[0][1]) / (4*ky-ky)

    mx = (4*kx-(kx-1))
    my = (4*ky - (ky-1))

    Setpoints = np.zeros(shape=(mx*my, 2))
    finite = np.zeros(shape=(kx*ky, 16))

    i2 = 0
    yi = 0
    num_points = 0
    num_elem = 0
    for j in range(my):
        xi = 0
        for i in range(mx):
            Setpoints[num_points][0] = x0+hx*i
            Setpoints[num_points][1] = y0+hy*j
            if j < my-1 and i < mx-1:
                if (xi % 3 == 0 and yi % 3 == 0):
                    tmp1 = [i2, i2 + 1, i2 + 2, i2 + 3, i2 + mx, i2 + mx + 1, i2 + mx + 2, i2 + mx + 3, i2 + mx + mx, i2 + mx + mx + 1, i2 + mx + mx + 2, i2 + mx + mx + 3, i2 + mx + mx + mx, i2 + mx + mx + mx + 1, i2 + mx + mx + mx + 2, i2 + mx + mx + mx + 3]
                    #tmp1 = [i2,i2+1,i2+mx,i2+mx+1]
                    finite[num_elem] = tmp1
                    num_elem = num_elem+1
            i2 = i2+1
            xi = xi+1
            num_points = num_points+1
        yi = yi+1


def splane(poin, femel):
    P = 0

    for i in range(len(femel)):
        x1 = Setpoints[femel[i][0]][0]
        x2 = Setpoints[femel[i][15]][0]
        hx = x2-x1
        y1 = Setpoints[femel[i][0]][1]
        y2 = Setpoints[femel[i][15]][1]
        hy = y2-y1
        if x1 <= poin[0] <= x2 and y1 <= poin[1] <= y2:
            num = i
            xsi = (poin[0]-x1)/hx
            nsi = (poin[1]-y1)/hy

    elem = femel[num]
    x1 = Setpoints[elem[0]][0]
    x2 = Setpoints[elem[15]][0]
    y1 = Setpoints[elem[0]][1]
    y2 = Setpoints[elem[15]][1]
    for i in range(16):
        P += x[elem[i]]*pfi_lag(i, poin[0],x1,x2,poin[1],y1,y2)
    return P


a = np.loadtxt('point.txt')

createmesh()
finite = np.int_(finite)

ig = [[] for i in range(len(Setpoints))]
ia = np.array([])
ja = np.array([])
gg = np.array([])
di = np.zeros(len(Setpoints))
x = np.zeros(len(Setpoints))
globalf = np.zeros(len(Setpoints))

create_ia_ja(len(finite), len(Setpoints))

num = 0

for elem in finite:
    b = np.zeros(shape=(16))
    matrix = np.zeros(shape=(16, 16))
    x1 = Setpoints[elem[0]][0]
    x2 = Setpoints[elem[15]][0]
    hx = x2 - x1
    y1 = Setpoints[elem[0]][1]
    y2 = Setpoints[elem[15]][1]
    hy = y2 - y1
    for i in range(16):
        for j in range(16):
            for k in range(len(a)):
                if x1 <= a[k][0] and a[k][0] <= x2 and y1 <= a[k][1] and a[k][1] <= y2:
                    matrix[i][j] += 1 * pfi_lag(i, a[k][0], x1, x2, a[k][1], y1, y2)*pfi_lag(j, a[k][0], x1, x2, a[k][1], y1, y2)

    for i in range(16):
        for k in range(len(a)):
            if x1 <= a[k][0] and a[k][0] <= x2 and y1 <= a[k][1] and a[k][1] <= y2:
                b[i] += 1 * pfi_lag(i, a[k][0], x1, x2, a[k][1], y1, y2)*a[k][2]

    A1 = np.zeros(shape=(16, 16))
    e = np.array([-1*sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7-2/7*sqrt(6/5)), -1*sqrt(3/7+2/7*sqrt(6/5)), sqrt(3/7+2/7*sqrt(6/5))])
    w = np.array([(18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36, (18-sqrt(30))/36])
    J = (x2-x1)*(y2-y1)/4
    for i in range(16):
        for j in range(16):
            for k in range(len(w)):
                for l in range(len(w)):
                    xsi = ((x2 + x1) / 2.0) + ((x2 - x1) / 2.0) * e[k]
                    nsi = ((y2 + y1) / 2.0) + ((y2 - y1) / 2.0) * e[l]
                    A1[i][j] += 1e-3 * w[k]*w[l]*d_pfi_lag_dx(i,xsi,x1,x2,nsi,y1,y2)*d_pfi_lag_dy(i,xsi,x1,x2,nsi,y1,y2)*d_pfi_lag_dx(j,xsi,x1,x2,nsi,y1,y2)*d_pfi_lag_dy(j,xsi,x1,x2,nsi,y1,y2)*J

    matrix = matrix+A1

    sborka(elem, matrix, b)

los_with_LU()

result_f = np.zeros(50)
result_x = np.zeros(50)
result_y = np.zeros(50)

tpoin = np.loadtxt('a.txt')

f = open('result.txt', 'w')
f1 = open('spline.txt', 'w')
i = 0
for item in a:
    P = splane(item, finite)
    norm = item[2]-P
    f1.write(str(item[0]) + " " + str(item[1]) + " "+ str(P) + '\n')
    for item_1 in tpoin:
        if (item[0] == item_1[0] and item[1] == item_1[1]):
            f.write(str(item[0]) + " " + str(item[1]) + " "+ str(item[2]) + " "+ str(P) + " "+ str(norm) + '\n')
    if item[1] == 1.569536423841062:
        result_f[i] = P
        result_x[i] = item[0]
        result_y[i] = item[2]
        i=i+1

print(i)
#X, Y = np.meshgrid(result_x, result_y)

# рисуем график
#plt.contour(X, Y, result_f)
#plt.plot(result_x, result_f)
fig = plt.figure(figsize=(6,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height])

plt.plot(result_x, result_y)
plt.plot(result_x, result_f)

ax.set_title('Сплайн, a=1e-3')
ax.set_xlabel('x')
ax.set_ylabel('P(x,y)')
# показываем график
plt.show()
