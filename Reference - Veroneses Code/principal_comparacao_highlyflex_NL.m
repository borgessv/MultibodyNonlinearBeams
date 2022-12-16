clear all;
addpath('multibody');
addpath('AeroFlex');

%% carrega parametros da asa
%wing = goland;
wing = highlyflex; 


%% inicializa e grava o modelo
nrb = 3; %numero de corpos rigidos
%matrizes_multibody = multibody(wing,nrb);
matrizes_multibody = multibody(wing,nrb,"coupled");
naf = 6; %numero de elementos aeroflex
matrizes_aeroflex = aeroflex(wing, naf, 0);

%% calculo aeroelastico

H = 20000; %altura das simulacoes
rho = atmosfera(H);
damp_ratio = 0.0001;
Dmb = matrizes_multibody.K*damp_ratio;
matrizes_aeroflex.CG = matrizes_aeroflex.KG*damp_ratio;


func =@(q) matrizes_multibody.K*q -matrizes_multibody.B_gravfun(q);% Equação de equilibrio
Xeq = fsolve(func,zeros(2*nrb,1)); 
X0multibody = [Xeq; zeros(2*nrb,1)];

[~, X0aeroflex] = trimairplane(matrizes_aeroflex,1,H,1,0,0);

V_vec = 3:0.5:10; 
%V_vec = 80:2:110
%V_vec = 6:0.5:10
i = 1;
mreal = 0;
Wn_mat = [];
Z_mat = [];
for V = V_vec
    Amultibody = aeroelasticity_multibody(X0multibody, V, rho, Dmb, matrizes_multibody);
    Aaeroflex = aeroelasticity_aeroflex(X0aeroflex',V,H, matrizes_aeroflex);
    %
    % Problema de Autovalor
    %[~, eigenValues] = eig(Q);
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
plot(NaN,NaN, 'bo');hold on;
plot(NaN,NaN,'ko');
plot(V_vec, Zmb_mat, 'bo');
plot(V_vec, Zaf_mat, 'ko');
legend('multibody', 'aeroflex')
grid on;
axis([min(V_vec) max(V_vec) -0.001 0.001]);
xlabel('velocidade');
ylabel('damping ratio');
subplot(212);
plot(V_vec, Wnmb_mat, 'bo'); hold on;
plot(V_vec, Wnaf_mat, 'ko');
xlabel('velocidade');
ylabel('frequencia (rad/s)');
axis([min(V_vec) max(V_vec) 0 40]);

%% plotar polos numa velocidade especifica
V = 7.2;
Amultibody = aeroelasticity_multibody(X0multibody, V, rho, Dmb, matrizes_multibody);
Aaeroflex = aeroelasticity_aeroflex(X0aeroflex',V,H, matrizes_aeroflex);

figure;
plot(NaN,NaN, 'bo');hold on;
plot(NaN,NaN,'ko');
plot(eig(Amultibody), 'bo'); 
plot(eig(Aaeroflex), 'ko'); 
legend('multibody', 'aeroflex')
axis([-2 0.5 -40 40]); grid on;
title('mapa de polos');
xlabel('real');
ylabel('imag');