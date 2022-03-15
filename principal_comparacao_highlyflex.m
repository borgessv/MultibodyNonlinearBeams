clear all;
addpath('galerkin');
addpath('multibody');
addpath('AeroFlex');

%% carrega parametros da asa
%wing = goland;
wing = highlyflex; 


%% inicializa e grava o modelo
matrizes_galerkin = galerkin(wing);
nrb = 3; %numero de corpos rigidos
%matrizes_multibody = multibody(wing,nrb);
matrizes_multibody = multibody(wing,nrb,"coupled");
naf = 3; %numero de elementos aeroflex
matrizes_aeroflex = aeroflex(wing, naf, 0);

%% calculo aeroelastico

H = 20000; %altura das simulacoes
rho = atmosfera(H);
damp_ratio = 0.0001;

Dgk = matrizes_galerkin.K*damp_ratio;
Dmb = matrizes_multibody.K*damp_ratio;
matrizes_aeroflex.CG = matrizes_aeroflex.KG*damp_ratio;

X0multibody = zeros(4*nrb,1);
X0aeroflex = zeros(length(matrizes_aeroflex.KG),1);

V_vec = 3:0.5:10; 
i = 1;
mreal = 0;
Wn_mat = [];
Z_mat = [];
for V = V_vec
    Agalerkin = aeroelasticity_galerkin(V, rho, Dgk, matrizes_galerkin);
    Amultibody = aeroelasticity_multibody(X0multibody, V, rho, Dmb, matrizes_multibody);
    Aaeroflex = aeroelasticity_aeroflex(X0aeroflex,V,H, matrizes_aeroflex);
    %
    % Problema de Autovalor
    %[~, eigenValues] = eig(Q);
    [Wn, Z] = damp(Agalerkin);
    Wn_mat(i,:) = Wn;
    Z_mat(i,:) = Z;
    [Wnmb, Zmb] = damp(Amultibody);
    Wnmb_mat(i,:) = Wnmb;
    Zmb_mat(i,:) = Zmb;
    [Wnaf, Zaf] = damp(Aaeroflex);
    Wnaf_mat(i,:) = Wnaf;
    Zaf_mat(i,:) = Zaf;
    %mreal(i) = max(real(diag(eigenValues)));
    i = i+1;
end
%%
figure;
subplot(211);
V_vec = V_vec;
plot(NaN,NaN, 'ro');hold on;
plot(NaN,NaN, 'bo');
plot(NaN,NaN,'ko');
plot(V_vec, Z_mat, 'ro'); 
plot(V_vec, Zmb_mat, 'bo');
plot(V_vec, Zaf_mat, 'ko');
legend('galerkin', 'multibody', 'aeroflex')
grid on;
axis([min(V_vec) max(V_vec) -0.001 0.001]);
xlabel('velocidade');
ylabel('damping ratio');
subplot(212);
plot(V_vec, Wn_mat, 'ro');hold on;
plot(V_vec, Wnmb_mat, 'bo');
plot(V_vec, Wnaf_mat, 'ko');
xlabel('velocidade');
ylabel('frequencia (rad/s)');
axis([min(V_vec) max(V_vec) 0 40]);

%% plotar polos numa velocidade especifica
V = 7.2;
Agalerkin = aeroelasticity_galerkin(V, rho, Dgk, matrizes_galerkin);
Amultibody = aeroelasticity_multibody(X0multibody, V, rho, Dmb, matrizes_multibody);
Aaeroflex = aeroelasticity_aeroflex(X0aeroflex,V,H, matrizes_aeroflex);

figure;
plot(NaN,NaN, 'ro');hold on;
plot(NaN,NaN, 'bo');
plot(NaN,NaN,'ko');
plot(eig(Agalerkin), 'ro'); hold on;
plot(eig(Amultibody), 'bo'); 
plot(eig(Aaeroflex), 'ko'); 
legend('galerkin', 'multibody', 'aeroflex')
axis([-2 0.5 -40 40]); grid on;
title('mapa de polos');
xlabel('real');
ylabel('imag');