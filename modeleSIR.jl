using LinearAlgebra
using PyPlot

Nx = 50
mid = convert(Int64,Nx/2)
Nt = 10
h = 1/Nt #pas de temps
D = 0.2#coefficient de diffusion
I = zeros(Nx,Nx,Nt)

#initialisation
I[mid, mid, 1] = 1

#schema d'euler explicite
Ay = Tridiagonal(ones(Nx-1), -2.0*ones(Nx), ones(Nx-1))/h^2
Ax = Tridiagonal(ones(Nx-1), -2.0*ones(Nx), ones(Nx-1))/h^2
Ay[1,2] =2.0
Ay[end, end-1] =2.0
Ax[2,1] = 2.0
Ax[end-1, end] =2.0
In = zeros(Nx,Nx)#variable de stockage
In[mid,mid] = 1

for t=1:Nt
    In[:,:] = In[:,:] + h*D*(Ay*In + In*Ax)
    I[:,:,t] = In[:,:]
end


#affichage
print(I[:,5,Nt])
figure()
X = [i for i=1:Nx]
Y = [i for i=1:Nx]
surf(X, Y, I[:,:,Nt], rstride=10, cstride=10)
plt.show()
plt.savefig()
