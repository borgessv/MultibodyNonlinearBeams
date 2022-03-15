function A = aeroelasticity_multibody(X0, U, rho, D, matrizes_multibody)
    A = lineariza(@(X)dinam_aerofun(X,matrizes_multibody,D,U, rho),X0); %Linearização
end 

function [result] = dinam_aerofun(X,est,D, U, rho)%,K_psi,I
n = est.n;
q = X(1:2*n,1);
p = X(2*n+1:4*n,1);
MM = est.MMfun(q);
qdot = MM\p;
pdot = -est.K*q - D*qdot + est.B_aerofun(q,qdot,U,rho);%+ B_grav(q);
result = [qdot;pdot];
end
