function [B] = gravity(q,w0,L0,a,d0)
    n = length(q)/2;
    Bdot = zeros(2*n,1);
    Bdot_T = zeros(2*n,1);
    
    for i = 1:n
        Bdot(i) = 1;
        %Bdot_T(n+i) = 1;
    end
    
    B = zeros(2*n,1);
    F = w0*L0;
    Ftip = w0*L0/2;
    M_z = w0*(L0^2)/12;
    nx = 0;
    ny = 0;
    nz = -1;
    for i = 1:n
        if i == n
            [~,~,~,~,~,jforca] = forceposition_3d(n,i,L0);
            Bmat = jforca*[nx;ny;nz];
            %Bdot_T = zeros(2*n,1);
            Bdot_T(n+i) = 1;
            B = B + Ftip*Bmat +F*(a-d0)*Bdot_T - M_z*Bdot;
        else
            [~,~,~,~,~,jforca] = forceposition_3d(n,i,L0);
            Bmat = jforca*[nx;ny;nz];
            %Bdot_T = zeros(2*n,1);
            Bdot_T(n+i) = 1;
            B = B + F*Bmat +F*(a-d0)*Bdot_T;
        end
    end
end