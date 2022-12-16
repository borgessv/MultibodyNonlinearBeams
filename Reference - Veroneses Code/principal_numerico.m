addpath("multibody_num");
wing = highlyflex;
n = 16;
xeq = fsolve(@(X) solve_equilibrium(X,wing), zeros(2*n,1));
%multibody_num(wing, n)

%%
h = elemental_position_num(xeq,16/n);
figure;plot3(h(1:n), h((n+1):(2*n)), -h((2*n+1):end))
axis equal
grid on
%%
matrizes = multibody_num(xeq, wing, n);

[eivec, eival] = (eig(matrizes.MM\matrizes.K));
eival = sqrt(diag(eival));
[eival ind] = sort(eival);
av_nl = eival(1:4);

matrizes = multibody_num(xeq*0, wing, n);
eival = sqrt(eig(matrizes.MM\matrizes.K));
eival = sort(eival);
av_l = eival(1:4);

[av_nl av_l]
