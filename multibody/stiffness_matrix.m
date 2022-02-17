function [K] = stiffness_matrix(n,L0,EI,GJ)
    disp = zeros(n,1);
    theta = zeros(n,1);
    dist = zeros(n,1);
    nload = n;
    P = 1;
    for i = 1:n
        disp(i) = P*(i^2)*(L0^3)*(3*n-i)/(6*EI);
    end

    dist(1) = disp(1);
    for i = 1:n-1
        dist(i+1) = disp(i+1) - disp(i);
    end

    theta(1) = dist(1)/L0;
    Mom = zeros(n,1);
    Mom(1) = P*n*L0;
    for i = 1:n-1
        theta(i+1) = (dist(i+1) - dist(i))/L0;
        Mom(i+1) = P*L0*(n-i) + 0*L0*(nload-i);
    end

    k_psi_n = GJ/(L0);
    k0_psi = ones(n,1)*k_psi_n;
    k0_psi(1) = k_psi_n*2;

    k0 = Mom./theta;
    K_theta = diag(k0);
    K_psi = diag(k0_psi);
    K = blkdiag(K_theta,K_psi);
end
