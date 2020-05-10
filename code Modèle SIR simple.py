import numpy as np
import matplotlib.pyplot as plt

## Paramétrages

# On prend les valeurs de beta, N et gamma valable pour le COVID-19 en France
beta = 1.03
N = 67*1e6
gamma = 0.0183

## définition de fonctions

def Euler_vect(F,X0,liste_t):
    N = len(liste_t)
    r = len(X0)
    X = np.zeros((N,r))
    X[0,:] = X0
    for i in range(N-1):
        h = liste_t[i+1] - liste_t[i]
        X[i+1,:] = X[i,:] + h*F(X[i,:],liste_t[i])
    return X

def f(X,t):
    S,I,R  = X[0],X[1],X[2]
    return np.array([-(beta/N)*I*S,(beta/N)*I*S -gamma*I, gamma*I])
    

    
## Tracés
# on étudie ce modèle sur 345 jours

#paramètres
N =1  #on normalise la population
Nb = 365 #nombre de jours d'étude
Pt = 1000 #nombre de points
t = np.linspace(0,Nb,Pt)

# resultats
X1 = Euler_vect(f, [N,1e-7,0],t)
X1_derivee = np.zeros((Pt,3))
for i in range(Pt):
    X1_derivee[i,:] = f(X1[t[i]],t[i])
S = [e[0] for e in X1] 
S_derivee = X1_derivee[:,0]
I = [e[1] for e in X1]
I_derivee = X1_derivee[:,1]
R = [e[2] for e in X1]
R_derivee = X1_derivee[:,2]

## tracés des solutions + portrait de phase
plt.plot(t,S)
plt.plot(t, I)
plt.plot(t,R)
## tracés portrait de phase
plt.plot(S,S_derivee,label = 'S')
plt.plot(I,I_derivee,label = 'I')
plt.plot(R,R_derivee,label = 'R')

plt.legend()

## tracé champ de vecteurs
N = 1 # on normalise la population en considérant I, S et R comme des fractions d'i-celle
[X,Y] = np.meshgrid (np.linspace(0,1,15) ,np.linspace(0,1,15))
DX = -(beta/N)*X*Y
DY= (beta/N)*X*Y - gamma*Y
V1 = DX/(DX**2+DY**2)**(1/2) 
V2=DY/(DX**2+DY**2)**(1/2)
plt.quiver(X,Y,V1,V2)

DX = -(beta/N)*X*Y
DY= (beta/N)*X*Y - gamma*Y
plt.contour(X,Y,DX)
plt.contour(X,Y,DY)
plt.show()
