function [L] = lift_qs(U,rho,c,a,L0,qq,qq_dot,z_dot)
alpha_0 = deg2rad(0);
CL = 2*pi*(qq_dot*(0.25*c+a)/U - z_dot/U +(qq - alpha_0));
L = 0.5*rho*c*CL*L0*U^2;
end