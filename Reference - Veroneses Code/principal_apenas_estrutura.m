clear all;
addpath('multibody');

%% carrega parametros da asa
wing = highlyflex; % Patil --> Hodges


%% inicializa e grava o modelo
nrb = 3; %numero de corpos rigidos
matrizes_multibody = multibody(wing,nrb,"coupled"); % a versao coupled
%acopla torcao e flexao e eh bem mais lenta para calcular a cinematica!
%matrizes_multibody = multibody(wing,nrb);

%%
% cálculo de equilíbrio
func =@(q) matrizes_multibody.K*q - matrizes_multibody.B_gravfun(q);% 
Xeq = fsolve(func,zeros(2*nrb,1)); 

h = matrizes_multibody.pos_(Xeq);
h = [0;h(1:nrb);0;h(nrb+1:2*nrb);0;h(2*nrb+1:3*nrb)];

figure;plot3(h(1:nrb+1), h((nrb+2):(2*nrb+2)), h((2*nrb+3):end)); axis equal;
xlabel('x')
ylabel('y')
zlabel('z')

x = h(1:nrb+1);
y = h((nrb+2):(2*nrb+2));
z = h((2*nrb+3):end);

for i = 1:length(x)
    v(:,4*i-3) = [x(i); y(i)+wing.b/2; z(i)+wing.b/2];
    v(:,4*i-2) = [x(i); y(i)+wing.b/2; z(i)-wing.b/2];
    v(:,4*i-1) = [x(i); y(i)-wing.b/2; z(i)-wing.b/2];
    v(:,4*i) = [x(i); y(i)-wing.b/2; z(i)+wing.b/2];
end
%%
%p = [p1'; p2'; p3'; p4'; p5'; p6'; p7'; p8'];
f = [1 5 9 13 14 10 6 2];

hold on
patch('Faces',f,'Vertices',v','FaceColor','b','edgecolor','none')