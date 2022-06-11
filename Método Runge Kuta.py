# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 18:45:29 2021
Método Numérico de Runge Kuta para el problema

@author: hiram
"""
from matplotlib import pyplot as plt

alpha = 0.0166
beta = 0.2875
lamda = 0.3895
mu = 0.0051
gamma = 0.1908
delta = 0.0008
rho = 0.3666
beta2 = 0.2625

No = 126
Do = 7.8144 
Co = 4.9856 
So = 32.6 
h = 0.01
t = 0

D = [Do]
C = [Co]
S = [So]
N = [No]
DUC = [Do+Co]
PRV = [(Do+Co)/No]
N2 = [No]

T = [0]

def fD(t, Di, Ci, Si, Ni):
    fD = ( 0.0166*0.2875 - 0.3895 - 0.0051)*Di+ ( 0.0166*0.2875 + 0.1908)*Ci + 0.2875*(Si*Di/Ni)
    return fD

def fC(t, Di, Ci):
    fC = 0.3895*Di - (0.1908 + 0.0008 + 0.0051)*Ci
    return fC

def fS(t, Di, Ci, Si, Ni):
    fS =  0.0166*Si +  0.0166*(1-0.3666)*Di +  0.0166*(1-0.3666)*Ci - 0.2875*(Si*Di/Ni)
    return fS

def fN(t,   Ni):
    fN = (0.0166 - 0.0051)*Ni
    return fN

def RK(t, Di, Ci, Si, Ni, h):
    "K1 -------------------------------------------------------------"
    K11 = fD(t, Di, Ci, Si, Ni)
    K21 = fC(t, Di, Ci)
    K31 = fS(t, Di, Ci, Si, Ni)
    K41 = fN(t, Ni)
    "K2 ------------------------------------------------------------"
    K12 = fD(t+h/2, Di + (K11*h)/2, Ci + (K21*h)/2, Si + (K31*h)/2, Ni + (K41*h)/2)
    K22 = fC(t+h/2, Di + (K11*h)/2, Ci + (K21*h)/2 )
    K32 = fS(t+h/2, Di + (K11*h)/2, Ci + (K21*h)/2, Si + (K31*h)/2, Ni + (K41*h)/2)  
    K42 = fN(t + h/2, Ni + (K41*h)/2) 
    "K3 --------------------------------------------------------------"
    K13 = fD(t+h/2, Di + (K12*h)/2, Ci + (K22*h)/2, Si + (K32*h)/2, Ni + (K42*h)/2)
    K23 = fC(t+h/2, Di + (K12*h)/2, Ci + (K22*h)/2 )
    K33 = fS(t+h/2, Di + (K12*h)/2, Ci + (K22*h)/2, Si + (K32*h)/2, Ni + (K42*h)/2)
    K43 = fN(t + h/2, Ni + (K42*h)/2)
    "K4  --------------------------------------------------------------"
    K14 = fD(t+h, Di + (K13*h), Ci + (K23*h), Si + (K33*h), Ni + (K43*h))
    K24 = fC(t+h, Di + (K13*h), Ci + (K23*h) )
    K34 = fS(t+h, Di + (K13*h), Ci + (K23*h), Si + (K33*h), Ni + (K43*h))
    K44 = fN(t + h, Ni + (K43*h))
    "Fi---------------------------------------------------------------"
    FD = Di + h*(K11+2*K12+2*K13+K14)/6
    FC = Ci + h*(K21+2*K22+2*K23+K24)/6
    FS = Si + h*(K31+2*K32+2*K33+K34)/6
    FN = Ni + h*(K41+2*K42+2*K43+K44)/6
    "------------------------------------------------------------------"
    return [FD, FC, FS, FN]

while t < 5:
    Kutta = RK(t, Do, Co, So, No, h)

    D.append(Kutta[0])
    C.append(Kutta[1])
    S.append(Kutta[2])
    N.append(Kutta[3])
    DUC.append(Kutta[0] + Kutta[1])
    PRV.append((Kutta[0] + Kutta[1])/Kutta[3])
    t += h
    T.append(t)
    Do = Kutta[0]
    Co = Kutta[1]
    So = Kutta[2]
    No = Kutta[3]


fig = plt.figure()
plt.plot(T, D, label= 'Diabéticos sin complicación')
plt.plot(T, C, label= 'Diabéticos con complicación')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title('Comparativa de diabéticos con y sin complicación')
plt.grid()
plt.legend()
plt.show()

fig2 = plt.figure()
plt.plot(T, DUC, label= 'Diabéticos')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title('Estimación de la pobalción diabética')
plt.grid()
plt.legend()
plt.show()

fig3 = plt.figure()
plt.plot(T, S, label= 'Suceptibles a la enfermedad')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title('Estimación de población suceptible')
plt.grid()
plt.legend()
plt.show()

fig4 = plt.figure()
plt.plot(T, N, label= 'Población')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title('Estimación de la población')
plt.grid()
plt.legend()
plt.show()

fig5 = plt.figure()
plt.plot(T, PRV, label= 'Prevalencia de la enfermedad')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title("Prevalencia de la diabetes en México")
plt.grid()
plt.legend()
plt.show()

print("Diabéticos sin complicaciones Inicial" , D[0], "Final" , D[-1], "Diferencia" , abs( D[-1]-D[0]))
print("Diabéticos con complicaciones Inicial" , C[0], "Final" , C[-1], "Diferencia" , C[-1]-C[0])
print("Suceptibles Inicial" , S[0], "Final" , S[-1], "Diferencia" , S[-1]-S[0])
print("Población Inicial" , N[0], "Final" , N[-1], "Diferencia" , N[-1]-N[0])
print("Diabéticos Inicial" , DUC[0], "Final" , DUC[-1], "Diferencia" , DUC[-1]-DUC[0])
print("Prevalencia Inicial" , PRV[0], "Final" , PRV[-1], "Diferencia" , PRV[-1]-PRV[0])