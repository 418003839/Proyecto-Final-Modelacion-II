# -*- coding: utf-8 -*-
"""
Método Numérico de Runge Kuta para el problema

@author: hiram
"""
from matplotlib import pyplot as plt
import math as mt

alpha = 0.0166
beta = 0.2875
lamda = 0.295
mu = 0.0051
gamma = 0.1908
delta = 0.0008
rho = 0.3666
beta2 = 0.2625

k=7
No = 126
Do = 7.8144 
Co = 4.9856 
So = 32.6 
h = 0.001
t = 0
t0=0

D = [Do]
C = [Co]
S = [So]
N = [No]
DUC = [Do+Co]
PRV = [(Do+Co)/No]


D1 = [Do]
C1 = [Co]
S1 = [So]
N1 = [No]
DUC1 = [Do+Co]
PRV1 = [(Do+Co)/No]


T = [0]
T1=[0]

def fD(t, Di, Ci, Si, Ni):
    fD = ( alpha*rho - lamda - mu)*Di+ ( alpha*rho + gamma)*Ci + beta*(Si*Di/Ni)
    return fD

def fC(t, Di, Ci):
    fC = lamda*Di - (gamma + delta+ mu)*Ci
    return fC

def fS(t, Di, Ci, Si, Ni):
    fS =  alpha*Si +  alpha*(1-rho)*(Di + Ci) - beta*(Si*Di/Ni)
    return fS

def fN(t,   Ni):
    fN = (alpha - mu)*Ni
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

while t<k:
    Kt = RK(t, Do, Co, So, No, h)

    D.append(Kt[0])
    C.append(Kt[1])
    S.append(Kt[2])
    N.append(Kt[3])
    DUC.append(Kt[0] + Kt[1])
    PRV.append((Kt[0] + Kt[1])/Kt[3])
    t += h
    T.append(t)
    Do = Kt[0]
    Co = Kt[1]
    So = Kt[2]
    No = Kt[3]
    if t-t0>=1:
        t0=mt.floor(t)
        D1.append(Kt[0])
        C1.append(Kt[1])
        S1.append(Kt[2])
        N1.append(Kt[3])
        DUC1.append(Kt[0] + Kt[1])
        PRV1.append((Kt[0] + Kt[1])/Kt[3])
        T1.append(t0)
    
        

fig = plt.figure()
plt.plot(T, D, label= 'Diabéticos sin complicación')
plt.plot(T, C, label= 'Diabéticos con complicación')
plt.scatter(T1, D1, label = 'Diabéticos sin complicación anual')
plt.scatter(T1, C1, label= 'Diabéticos con complicación anual')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title('Comparativa de diabéticos con y sin complicación')
plt.grid()
plt.legend()
plt.show()

fig2 = plt.figure()
plt.plot(T, DUC, label= 'Diabéticos')
plt.scatter(T1, DUC1, label= 'Diabéticos anuales')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title('Estimación de la pobalción diabética')
plt.grid()
plt.legend()
plt.show()

fig3 = plt.figure()
plt.plot(T, S, label= 'Suceptibles a la enfermedad')
plt.scatter(T1, S1, label= 'Suceptibles a la enfermedad anuales')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title('Estimación de población suceptible')
plt.grid()
plt.legend()
plt.show()

fig4 = plt.figure()
plt.plot(T, N, label= 'Población')
plt.scatter(T1, N1, label= 'Población anual')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title('Estimación de la población')
plt.grid()
plt.legend()
plt.show()

fig5 = plt.figure()
plt.plot(T, PRV, label= 'Prevalencia de la enfermedad')
plt.scatter(T1, PRV1, label= 'Prevalencia de la enfermedad anual')
plt.xlabel("Tiempo (Años)")
plt.ylabel("Población en millones de personas")
plt.title("Prevalencia de la diabetes en México")
plt.grid()
plt.legend()
plt.show()


"""
Esta pequeña parte del código genera una lista en excel donde da los datos anuales
para poder tenerlos a la mano. El problema de la creación del archivo radica en la 
ubicación de la creacion, así que si se desea el archivo, solo reemplazar 
después del excel_writer poner "C:/Users/"nombre de usuario"/"el resto de la ubicación"

Puedes tambien en la carpeta copiar su dirección y pegarla, sambiando la diagonal 
inversa (\) por una diagonal normal(/)
"""
"""
import pandas as pd

array = [T1,D1,C1,S1,N1,DUC1,PRV1]

df = pd.DataFrame(array).T
df.to_excel(excel_writer = "C:/Users/hiram/Desktop/Expo Modelación/Martínez Aragón Hiram/Diab2.xlsx")
"""
