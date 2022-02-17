function [result] = dinam_aero_red(t,Xr,nred,Kn,U)%,K_psi,I
q = Xr(1:nred,1);
p = Xr(nred+1:2*nred,1);
qdot = MMfunred_tor(q)\p;
C = 1e-4*Kn;
pdot = -Kn*q - C*qdot +B_aero_red(q,qdot,U);% + B_grav_red(q);
result = [qdot;pdot];
end