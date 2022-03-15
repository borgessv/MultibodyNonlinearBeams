%% STATIC SOLUTION

clc;clear;%close all

% Symbolics variables creation---------------------------------------------
n = 6;
nred = 2;
syms q [2*n 1]
syms q_dot [2*n 1]
syms p [2*n 1]
syms x [n 1];
syms y [n 1];
syms z [n 1];
syms U;

% Goland Wing--------------------------------------------------------------
L0 = 6.096/n; %Comprimento do elemento
c = 1.8288; %Corda (Largura do elemento)
a = 0.17*c; %Distância do eixo elastico em relação a metade da corda
m0 = 35.7186*L0; %Massa do elemento
rho = 1.225; %Densidade do ar
d0 = 0.06*c; %Distância do centro de massa em relação a metade da corda
I0 = 8.64*L0; %Inércia de cada elemento
EI = 9.77e6; %Rigidez de flexão
GJ = 0.987e6; %Rigidez de torção

% Peters very flexible wing-------------------------------------------------
% L0 = 16/n; 
% c = 1;
% a = 0;
% m0 = 0.75*L0;
% rho = 0.0889;
% d0 = 0;
% I0 = 0.1*L0; 
% F = 0;
% EI = 2e4;
% GJ = 1e4;

% Global coordinates calculation-------------------------------------------
% Calcula as posições dos elementos (centro de massa) nas coordenadas
% globais x,y,z em função das coordenadas generalizadas q
[pos,Jpos] = elemental_position(n,L0,d0,a);

%Stiffness matrix----------------------------------------------------------
%Calcula as rigidezes das molas com base na rigidez da viga
[K] = stiffness_matrix(n,L0,EI,GJ);

% System matrices----------------------------------------------------------

M = m0*eye(3*n); % Matriz de massa linear
I = I0*eye(2*n); % Matriz de inércia
T = zeros(n); %
T_1 = tril(ones(n)); %
T = blkdiag(T,T_1); % Matriz de transposição para incluir a matriz de inércia
X0 = zeros(2*n,1); % Deflexões iniciais
MM = transpose(Jpos)*M*Jpos + transpose(T)*I*T; %Cálculo da matriz de massa
%não linear
MMvec = double(subs(MM,q,X0)); %Matriz de massa não linear em repouso

% Natural modes------------------------------------------------------------
[eivec, eival] = eig(MMvec\K); %Autovalores e autovetores
[d,ind] = sort(diag(eival));
eigenValues = sqrt(diag(eival(ind,ind)));
eigenVectors = eivec(:,ind);

%Distributed force---------------------------------------------------------
U = 0; %Velocidade para simulação estática
w0 = 9.81*m0/L0; % Força distribuida de gravidade
B_grav = gravity(q,w0,L0,a,d0); %Calculo do vetor não linear de forças
%generalizadas geradas pela gravidade
[B_aero,L_vec] = aero_loads(rho,U,L0,q,q_dot,c,a); %Calculo do vetor
%não linear de forças generalizadas geradas pela aerodinamica
B_aero = subs(B_aero,q_dot,zeros(length(q_dot),1)); % Zerando as velocidades
%para fazer a simulação estática

% Reduced model--------------------------
syms eta [nred 1]
syms eta_dot [nred 1]
V = eigenVectors(:,1:nred); %Autovetores do modelo reduzido
Vi = V*eta;
Vi_dot = V*eta_dot;
Kn = transpose(V)*K*V; %Matriz de rigidez reduzida
pp = subs(pos,q,Vi); %posição globais em termo da variavel eta (Modelo
%reduzido)
Jpeta = jacobian(pp,eta);
MMred = transpose(Jpeta)*M*Jpeta; %Matriz de massa não linear reduzida 
B_grav_red = transpose(V)*subs(B_grav,q,Vi);
B_aero_red = transpose(V)*subs(B_aero,q,Vi);
% Solution calculation-----------------------------------------------------
tic()
f = K*q -B_grav - B_aero;% Equação de equilibrio
func = matlabFunction(f,'Vars',{q});
Xeq = fsolve(func,X0); 
toc()

% Solution calculation red-------------------------------------------------
eta0 = zeros(nred,1);
tic()
fred = Kn*eta -B_grav_red - B_aero_red;%
funcred = matlabFunction(fred,'Vars',{eta});
etaeq = fsolve(funcred,eta0);
toc()

L_plot = subs(L_vec,q,Xeq); %Esse Lplot ele recupera da função aerodinamica
%para conseguir plotar
L_plot = subs(L_plot,q_dot,Xeq*0);
plot_struct_3d_aero(L0,c,a,L_plot,Xeq); %Essa função plota a estrutura
hold on
plot_struct_3d_aero(L0,c,a,L_plot,V*etaeq); %Essa função plota a estrutura
% (modelo reduzido)
Lift = double(sum(L_plot)); %Verificação da sustentação total
%% DYNAMIC SOLUTION
% Functions creation:
% A equação dinamica está dentro da função dinam_aero e dinam_aero_red para
% o modelo reduzido, como essas funções para o calculo da matriz de massa
% não linear e aerodinamica são com matemática simbólica, para usar na
% simulação é bem mais eficiente criar essas funções como arquivos usando o
% matlabFunction. Então toda a vez que mexer no modelo da asa ou nas
% funções que geram essas funções, é necessário rodar essa bloco para criar
% essas funções novamente. Como mexemos bastante no valor da velocidade,
% para não ter que criar essas funções toda a vez que mudar, a velocidade é
% colocada como input dessas funções
syms U
[B_aero,L_vec] = aero_loads(rho,U,L0,q,q_dot,c,a); %Calculo do vetor
%não linear de forças generalizadas geradas pela aerodinamica com
%velocidade de input
functio_red = subs(B_aero,q,Vi);
functio_red2 = subs(functio_red,q_dot,Vi_dot); 
B_aero_red = transpose(V)*functio_red2; % Função aerodinamica para o modelo
%reduzido
matlabFunction(B_aero,'File','B_aero','Vars',{q;q_dot;U})
matlabFunction(B_grav,'File','B_grav','Vars',{q})
%matlabFunction(MM,'File','MMfun_tor','Vars', {q})
%matlabFunction(B_aero_red,'File','B_aero_red','Vars',{eta;eta_dot;U})
%matlabFunction(B_grav_red,'File','B_grav_red','Vars',{eta})
%matlabFunction(MMred,'File','MMfunred_tor','Vars', {eta})
%% AEROELASTIC ANALYSIS
i = 1;
U_vec = (80:2:110)*0.3048;
for U = U_vec
    X0 = [Xeq*0; Xeq*0]; %Estado inicial
    t = 0;
    A = lineariza(@(X)dinam_aero(t,X,n,K,U),X0); %Linearização

    [Wn, Z] = damp(A);
    Wn_mat(i,:) = Wn;
    Z_mat(i,:) = Z;
    %mreal(i) = max(real(diag(eigenValues)));
    i = i+1;
end

figure;
subplot(211);
V_vec = U_vec;
plot(V_vec, Z_mat, 'ro'); hold on;
axis([min(V_vec) max(V_vec) -0.001 0.001]);
subplot(212);
plot(V_vec, Wn_mat, 'ro'); hold on;
axis([min(V_vec) max(V_vec) 0 100]);
%U = 25; %Velocidade
%%
U = 20;
figure;
stability_plot(A)
%plot de estabilidade

%% DYNAMIC SOLUTION REDUCED
U = 10;
eta0 = [etaeq*0;etaeq*0];
t = 0;
Ar = lineariza(@(Xr)dinam_aero_red(t,Xr,nred,Kn,U),eta0);
stability_plot(Ar)