clear all;
addpath('multibody');

%% carrega parametros da asa
wing = highlyflex; % Patil --> Hodges


%% inicializa e grava o modelo
nrb = 3; %numero de corpos rigidos
%matrizes_multibody = multibody(wing,nrb,"coupled"); % a versao coupled
%acopla torcao e flexao e eh bem mais lenta para calcular a cinematica!
matrizes_multibody = multibody(wing,nrb);

%%
% cálculo de equilíbrio
func =@(q) matrizes_multibody.K*q - matrizes_multibody.B_gravfun(q);% 
Xeq = fsolve(func,zeros(2*nrb,1)); 

h = matrizes_multibody.pos_(Xeq);
h = [0;h(1:nrb);0;h(nrb+1:2*nrb);0;h(2*nrb+1:3*nrb)];

figure;plot3(h(1:nrb+1), h((nrb+2):(2*nrb+2)), h((2*nrb+3):end)); axis equal;

!git add .
!git commit -m "add existing file"
!git push origin main
