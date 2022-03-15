function [pos,Jpos,J_z] = position_elastic(n,r,L0)
    syms q [2*n 1]
    qq = q(1);
    x_b = 0.5*L0*cos(qq);
    y_b = 0;
    z_b = 0.5*L0*sin(qq);
    for i = 2:r
        x_b = x_b + 0.5*L0*cos(qq) + 0.5*L0*cos(qq + q(i));
        y_b = 0;
        z_b = z_b + 0.5*L0*sin(qq) + 0.5*L0*sin(qq + q(i));
        qq = qq + q(i);
    end
    x = x_b;
    y = y_b;
    z = z_b;
    pos = [x;y;z];
    Jpos = jacobian(pos,q);
    J_z = jacobian(z,q);
end