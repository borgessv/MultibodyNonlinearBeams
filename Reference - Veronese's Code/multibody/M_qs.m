function [M_theta] = M_qs(U,rho,L0,c,qq_dot)
    M_theta = qq_dot*L0*(pi*rho*U*(c/2)^3)/2;
end