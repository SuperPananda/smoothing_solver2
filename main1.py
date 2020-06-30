import numpy as np
from scipy import sparse
from math import sqrt, fabs, log
import matplotlib.pyplot as plt

def pfi(i,x,x1,x2):
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


A = np.loadtxt('spline.txt')
b = np.zeros(shape=(len(A), 3))

""" j=0
for i in range(len(A)):
    if 1<=A[i][0] and A[i][0]<=2 and 1<=A[i][1] and A[i][1]<=2:
        b[j][0] = A[i][0]
        b[j][1] = A[i][1]
        b[j][2] = A[i][2]
        j=j+1

del_arr = np.delete(b, np.where(b == [0, 0, 0]), axis=0) """
#print(del_arr)

""" num = 0
for item in del_arr:
    if item[1] == 1.1125827814569538:
        num = num+1 """

""" f = open('result_a.txt', 'w')
for item in del_arr:
    f.write(str(item[0]) + " " + str(item[1]) + " "+ str(item[2]) + '\n') """

x_value = np.array(A[47:97, 0])
y_value = np.array(A[47:96, 0])

X, Y = np.meshgrid(x_value, y_value)

Z = np.zeros(shape=(len(Y),len(X[0])))

for k in range(len(A)):
    for i in range(len(Y)):
        for j in range(len(Y[0])):
            if (A[k][1] == Y[i][0]) and (A[k][0]==X[i][j]):
                Z[i][j] = A[k][2]       

fig = plt.figure(figsize=(6,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height]) 

cp = plt.contourf(X, Y, Z, cmap=plt.cm.RdYlGn)
plt.colorbar(cp)

ax.set_title('Сглаженное векторное решение, a=1e-3')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()



print(b)

node = np.array([0, 0.5, 1.1, 1,5, 2])

finite = np.array([0, 2])

e = np.array([-1*sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7-2/7*sqrt(6/5)), -1*sqrt(3/7+2/7*sqrt(6/5)), sqrt(3/7+2/7*sqrt(6/5))])
w = np.array([(18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36, (18-sqrt(30))/36])
h = 2
A = np.zeros(shape=(4, 4))
A1 = np.array([[2, -2, 1, -1],[-2, 2, -1, 1],[1, -1, 2, -2],[-1, 1, -2, 2]])
A2 = np.array([[2, 1, -2, -1],[1, 2, -1, -2],[-2, -1, 2, 1],[-1, -2, 1, 2]])
A3 = 2/2*1/6*A1+2/2*1/6*A2

x1 = 0
x2 = 2
y1 = 0
y2 = 2
J = (x2-x1)*(y2-y1)/4
for i in range(4):
    for j in range(4):
        for k in range(len(w)):
            for l in range(len(w)):
                xsi = ((x2 + x1) / 2.0) + ((x2 - x1) / 2.0) * e[k]
                nsi = ((y2 + y1) / 2.0) + ((y2 - y1) / 2.0) * e[l]
                A[i][j] += 1 * w[k]*w[l]*(d_pfi_dx(i,xsi,x1,x2,nsi,y1,y2)*d_pfi_dx(j,xsi,x1,x2,nsi,y1,y2)+d_pfi_dy(i,xsi,x1,x2,nsi,y1,y2)*d_pfi_dy(j,xsi,x1,x2,nsi,y1,y2))*J

print(A3)
print(A)
print(A)
