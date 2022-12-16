function [pos, Jpos] = elemental_position(n,L0,d0,a,varargin)
    if nargin > 4
        flag = varargin{1};
    else
        flag = "false";
    end
    syms q [2*n 1]
    syms x [n 1];syms x_T [n 1]
    syms y [n 1];syms y_T [n 1]
    syms z [n 1];syms z_T [n 1]
    
    if flag == "coupled"
    %     qq = q(1);
        rx = [1;0;0];
        ry = [0;-1;0];
        rz = [0;0;1];

        %rotacao em torno de x, angulo de torcao q(n+1)
        rx = rx;
        ry = rotationmatrix(rx, q(n+1))*ry;
        rz = rotationmatrix(rx, q(n+1))*rz;
        %rotacao em torno de y, angulo de flexao q(1)
        rx = rotationmatrix(ry, q(1))*rx;
        ry = ry;
        rz = rotationmatrix(ry, q(1))*rz;
        rx = simplify(rx);
        ry = simplify(ry);
        rz = simplify(rz);
        x(1) = rx(1)*0.5*L0;
        y(1) = rx(2)*0.5*L0;
        z(1) = rx(3)*0.5*L0;

        for i = 2:n
            x(i) = x(i-1) + 0.5*L0*rx(1);
            y(i) = y(i-1) + 0.5*L0*rx(2);
            z(i) = z(i-1) + 0.5*L0*rx(3);
            %rotacao em torno de x, angulo de torcao q(n+1)
            rx = rx;
            ry = rotationmatrix(rx, q(n+i))*ry;
            rz = rotationmatrix(rx, q(n+i))*rz;
            %rotacao em torno de y, angulo de flexao q(1)
            rx = rotationmatrix(ry, q(i))*rx;
            ry = ry;
            rz = rotationmatrix(ry, q(i))*rz;
            %rx = simplify(rx);
            %ry = simplify(ry);
            %rz = simplify(rz);
            x(i) = x(i) + 0.5*L0*rx(1);
            y(i) = y(i) + 0.5*L0*rx(2);
            z(i) = z(i) + 0.5*L0*rx(3);
        end
        pos = simplify([x;y;z]);        
    else
        
         qq = q(1);
         x(1) = 0.5*L0*cos(qq);
         y(1) = 0;
         z(1) = 0.5*L0*sin(qq);

         for i = 2:n
             x(i) = x(i-1) + 0.5*L0*cos(qq) + 0.5*L0*cos(qq + q(i));
             y(i) = 0;
             z(i) = z(i-1) + 0.5*L0*sin(qq) + 0.5*L0*sin(qq + q(i));        
             qq = qq + q(i);
         end
         qq = 0;
        qq_T = 0;
        for i = 1:n
            qq = qq + q(i);
            qq_T = qq_T + q(n+i);
            x_T(i) = -(d0-a)*sin(qq)*sin(qq_T);
            y_T(i) = (d0-a)*cos(qq_T);
            z_T(i) = (d0-a)*sin(qq_T)*cos(qq);
        end
        pos = [x+x_T;y+y_T;z+z_T];

    end

    Jpos = jacobian(pos,q);
end

function R = rotationmatrixs(vec, theta)
    ux = vec(1); uy = vec(2); uz = vec(3);
    R = [1, -uz*(theta), +uy*(theta);
        +uz*(theta), 1, -ux*(theta);
        -uy*(theta), +ux*(theta), 1];
end

function R = rotationmatrix(vec, theta)
    ux = vec(1); uy = vec(2); uz = vec(3);
    R = [cos(theta)+ux^2*(1-cos(theta)), ux*uy*(1-cos(theta))-uz*sin(theta), ux*uz*(1-cos(theta))+uy*sin(theta);
        uy*ux*(1-cos(theta))+uz*sin(theta), cos(theta)+uy^2*(1-cos(theta)), uy*uz*(1-cos(theta))-ux*sin(theta);
        uz*ux*(1-cos(theta))-uy*sin(theta), uz*uy*(1-cos(theta))+ux*sin(theta), cos(theta)+uz^2*(1-cos(theta))];
end