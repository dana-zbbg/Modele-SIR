using LinearAlgebra
using PyPlot

Nx = 50
mid = convert(Int64,Nx/2)
Nt = 10
h = 1/Nt #pas de temps
D = 1#coefficient de diffusion
I = zeros(Nx,Nt)

#initialisation
I[mid, 1] = 0.1

#schema d'euler implicite
A = Tridiagonal(ones(Nx-1), -2.0*ones(Nx), ones(Nx-1))/h^2
Id = Diagonal(ones(Nx))
In = zeros(Nx)#variable de stockage
In[mid] = 1
B = inv(Id-h*D*A)
for t=1:Nt
    In[:] = B*In
    I[:,t] = In[:]
end


#affichage
print(I[:,Nt])
figure()
X = [i for i=1:Nx]
T = [t for t=1:Nt]
surf(T, X, I, rstride=10, cstride=10)
plt.show()
plt.savefig()
