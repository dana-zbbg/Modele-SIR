using LinearAlgebra
using PyPlot

Nx = 100
dx = 1/Nx
mid = convert(Int64,Nx/2)
h = 1e-3 #pas de temps
Nt = 100
D = 1#coefficient de diffusion
gamma = 1 #coefficient de guérison
beta =0.1 #coefficient de contamination
S = zeros(Nx, Nt)
I = zeros(Nx,Nt)
R = zeros(Nx, Nt)

#initialisation
I[mid, 1] = 0.2
I[mid-1, 1] = 0.2
I[mid+1, 1] = 0.2
S[:,1] = ones(Nx)
S[mid, 1] = 0.8
S[mid+1, 1] = 0.8
S[mid-1, 1] = 0.8

#construction des matrices et des vecteurs du schema
A = Tridiagonal(ones(3*Nx-1), -2.0*ones(3*Nx), ones(3*Nx-1))*Nx^2
A[Nx+1,Nx] =0
A[Nx,Nx+1] = 0
A[2*Nx+1,2*Nx] =0
A[2*Nx,2*Nx+1] = 0
Id = Diagonal(ones(3*Nx))
B = zeros(3*Nx, 3*Nx)
for i=Nx+1:2*Nx
    B[i,i] = 1
    B[i+Nx,i] = -1
end

#initialisation variables de stockage
SIR = zeros(3*Nx)#stockage des grangeurs S,I et R
SI = zeros(3*Nx) #produit SxI
SIR[1:Nx] = S[:,1]
SIR[Nx+1:2*Nx] = I[:,1]
SIR[2*Nx+1:3*Nx] = R[:,1]
for i=1:Nx
    SI[i] = S[i,1]*R[i,1]
    SI[Nx+i] = -SI[i]
end
C = inv(Id-h*D*A+h*gamma*B)

#implementation du schema Euler implicite
for t=2:Nt
    SIR[:] = C*(SIR[:]-h*beta*SI[:])
    S[:,t] = SIR[1:Nx]
    I[:,t] = SIR[Nx+1:2*Nx]
    R[:,t] = SIR[2*Nx+1:3*Nx]
    for i=1:Nx
        SI[i] = S[i,t]*R[i,t]
        SI[Nx+i] = -SI[i]
    end
end


#affichage
figure()
X = [i*dx for i=1:Nx]
T = [t*h for t=1:Nt]
surf(T, X, S, rstride=10, cstride=10)
surf(T, X, I, rstride=10, cstride=10)
surf(T, X, R, rstride=10, cstride=10)
xlabel("temps")
ylabel("espace")
zlabel("densité de population")
plt.show()
#plt.savefig()
