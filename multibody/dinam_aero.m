function [result] = dinam_aero(t,X,n,K,U)%,K_psi,I
q = X(1:2*n,1);
p = X(2*n+1:4*n,1);
MM = MMfun_tor(q);
qdot = MM\p;
C = 1e-3*K + 1e-3*MM;
pdot = -K*q - C*0*qdot + B_aero(q,qdot,U);%+ B_grav(q);
result = [qdot;pdot];
end
