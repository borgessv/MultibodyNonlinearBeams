function saida = solve_equilibrium(X, wing)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
matrizes = multibody_num(X, wing, length(X)/2);
saida = matrizes.K * X - matrizes.grav;
end

