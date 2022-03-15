function [pos,Jpos,J_z] = aeroposition_3d(n,r,L0,c,a)
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
    qq = 0;
    qq_T = 0;
    for i = 1:r
        qq = qq + q(i);
        qq_T = qq_T + q(n+i);
        x_T = -(0.25*c-a)*sin(qq_T)*sin(qq);
        y_T = (0.25*c-a)*cos(qq_T);
        z_T = (0.25*c-a)*sin(qq_T)*cos(qq);
    end
    x = x_b+ x_T;
    y = y_b+ y_T;
    z = z_b+ z_T;
    pos = [x;y;z];
    Jpos = jacobian(pos,q);
    J_z = jacobian(z,q);
end