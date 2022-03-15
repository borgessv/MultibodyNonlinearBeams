function matrizes = multibody_num(X, wing, n)

L0 = wing.L/n; 
c = wing.b*2;
a = -wing.a*c/2;
m0 = wing.m*L0;
d0 = -wing.xa*c/2 + a;

I0 = wing.I*L0; 
EI = wing.EI;
GJ = wing.GJ;

% Global coordinates calculation-------------------------------------------
% Calcula as posições dos elementos (centro de massa) nas coordenadas
% globais x,y,z em função das coordenadas generalizadas q
%[pos,Jpos] = elemental_position(n,L0,d0,a,flag);

%Stiffness matrix----------------------------------------------------------
%Calcula as rigidezes das molas com base na rigidez da viga
matrizes.K = stiffness_matrix(n,L0,EI,GJ);

% System matrices----------------------------------------------------------

M = m0*eye(3*n); % Matriz de massa linear
I = I0*eye(2*n); % Matriz de inércia
T = zeros(n); %
T_1 = tril(ones(n)); %
T = blkdiag(T,T_1); % Matriz de transposição para incluir a matriz de inércia
X0 = zeros(2*n,1); % Deflexões iniciais
Jpos = lineariza_complexstep(@(q) elemental_position_num(q,L0), X);
matrizes.MM = transpose(Jpos)*M*Jpos + transpose(T)*I*T; %Cálculo da matriz de massa

%eival = sqrt(eig(MM\K));
%eival = sort(eival);
%eival(1:4)
w0 = 9.81*m0; % Força distribuida de gravidade

matrizes.grav = Jpos((2*n+1):end,:)'*ones(n,1)*w0;
end

