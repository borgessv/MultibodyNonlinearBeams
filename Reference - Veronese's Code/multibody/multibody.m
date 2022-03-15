function matrizes = multibody(wing, n, varargin)
    if nargin > 2
        flag = varargin{1};
    else
        flag = "false";
    end

syms q [2*n 1]
syms q_dot [2*n 1]
syms p [2*n 1]
syms x [n 1];
syms y [n 1];
syms z [n 1];
syms U;


L0 = wing.L/n; 
c = wing.b*2;
a = -wing.a*c/2;
m0 = wing.m*L0;
d0 = -wing.xa*c/2 + a;

I0 = wing.I*L0; 
EI = wing.EI;
GJ = wing.GJ;
% L0 = 6.096/n; %Comprimento do elemento
% c = 1.8288; %Corda (Largura do elemento)
% a = 0.17*c; %Distância do eixo elastico em relação a metade da corda
% m0 = 35.7186*L0; %Massa do elemento
% rho = 1.225; %Densidade do ar
%d0 = 0.06*c; %Distância do centro de massa em relação a metade da corda
% I0 = 8.64*L0; %Inércia de cada elemento
%I0 = I0*1.1572;
% EI = 9.77e6; %Rigidez de flexão
% GJ = 0.987e6; %Rigidez de torção

% Global coordinates calculation-------------------------------------------
% Calcula as posições dos elementos (centro de massa) nas coordenadas
% globais x,y,z em função das coordenadas generalizadas q
[pos,Jpos] = elemental_position(n,L0,d0,a,flag);

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

w0 = 9.81*m0/L0; % Força distribuida de gravidade
B_grav = gravity(q,w0,L0,a,d0); %Calculo do vetor não linear de forças
%generalizadas geradas pela gravidade
syms U rho
[B_aero,L_vec] = aero_loads(rho,U,L0,q,q_dot,c,a); %Calculo do vetor

posfun = matlabFunction(pos, 'Vars',{q});
B_aerofun = matlabFunction(B_aero,'Vars',{q;q_dot;U; rho});
B_gravfun = matlabFunction(B_grav,'Vars',{q});
MMfun = matlabFunction(MM,'Vars', {q});

matrizes.pos_ = posfun;
matrizes.K = K;
matrizes.M_ = MM;
matrizes.B_aerofun = B_aerofun;
matrizes.B_gravfun = B_gravfun;
matrizes.MMfun = MMfun;
matrizes.n = n;
end

