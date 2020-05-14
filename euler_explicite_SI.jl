using LinearAlgebra
using PyPlot

Nx = 100
dx = 1/Nx
mid = convert(Int64,Nx/2)
h = 1e-5 #pas de temps
Nt = 40000
D = 0#coefficient de diffusion
gamma = 0.01 #coefficient de guérison
beta = 10 #coefficient de contamination
S = zeros(Nx, Nt)
I = zeros(Nx,Nt)

#initialisation
for x=1:Nx
    S[x,1] = sin(x*pi*1/Nx)
    #I[x,1] = 0.5*sin(x*pi*1/Nx)
end
I[mid,1] = 0.2
I[mid-1,1] = 0.1
I[mid+1,1] = 0.1
I[mid-2,1] = 0.05
I[mid+2,1] = 0.05
S[:,1] = S[:,1]-I[:,1]

#construction des matrices et des vecteurs du schema
A = Tridiagonal(ones(Nx-1), -2.0*ones(Nx), ones(Nx-1))*Nx^2

#initialisation variables de stockage
SI = zeros(Nx) #produit SxI
for i=1:Nx
    SI[i] = S[i,1]*I[i,1]
end

#implementation du schema Euler implicite
for t=2:Nt
    S[:,t] = S[:,t-1]+h*D*A*S[:,t-1]-h*beta*SI[:]
    I[:,t] = I[:,t-1]+h*D*A*I[:,t-1]+h*beta*SI[:]-h*gamma*I[:,t-1]
    for i=1:Nx
        SI[i] = S[i,t]*I[i,t]
    end
end


#affichage
figure()
X = [i*dx for i=1:Nx]
T = [t*h for t=1:Nt]
surf(T, X, S, rstride=10, cstride=10)
surf(T, X, I, rstride=10, cstride=10)
xlabel("temps")
ylabel("espace")
zlabel("densité de population")
plt.show()
#plt.savefig()
