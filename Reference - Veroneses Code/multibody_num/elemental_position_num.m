function [pos] = elemental_position_num(q,L0)

	n = length(q)/2;
    x = zeros(n,1);
    y = zeros(n,1);
    z = zeros(n,1);
    
    rx = [1;0;0];
    ry = [0;1;0];
    rz = [0;0;1];

    %rotacao em torno de x, angulo de torcao q(n+1)
    rx = rx;
    ry = rotationmatrix(rx, q(n+1))*ry;
    rz = rotationmatrix(rx, q(n+1))*rz;
    %rotacao em torno de y, angulo de flexao q(1)
    rx = rotationmatrix(ry, q(1))*rx;
    ry = ry;
    rz = rotationmatrix(ry, q(1))*rz;
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
    pos = ([x;y;z]);        
    
   
end


function R = rotationmatrix(vec, theta)
    ux = vec(1); uy = vec(2); uz = vec(3);
    R = [cos(theta)+ux^2*(1-cos(theta)), ux*uy*(1-cos(theta))-uz*sin(theta), ux*uz*(1-cos(theta))+uy*sin(theta);
        uy*ux*(1-cos(theta))+uz*sin(theta), cos(theta)+uy^2*(1-cos(theta)), uy*uz*(1-cos(theta))-ux*sin(theta);
        uz*ux*(1-cos(theta))-uy*sin(theta), uz*uy*(1-cos(theta))+ux*sin(theta), cos(theta)+uz^2*(1-cos(theta))];
end