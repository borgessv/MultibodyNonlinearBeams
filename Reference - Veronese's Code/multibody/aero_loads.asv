function [B,L_vec] = aero_loads(rho,U,L0,q,q_dot,c,a)
n = length(q)/2; %Numero de elementos
Bdot_T = zeros(2*n,1); %Vetor unitario que contabiliza o momento de torção,
%no calculo análitico é o jacobiano para as torções
B = zeros(2*n,1); %Vetor de saída contendo as forças generalizadas
qq_Tdot = 0; %Angulo acumulado de torção
qq_T = 0; %Angulo acumulado de flexão
L_vec = sym(zeros(n,1)); %Vetor de sustentação
qq_dot = sym(zeros(2*n,1)); %Velocidade angular acumulada de flexao
for i = 1:n
    qq_T = qq_T+q(n+i); %Adiciona o angulo de torção do elemento
    qq_Tdot = qq_Tdot+q_dot(n+i);%Adiciona a velocidade angular de torção
    %do elemento
    qq_dot(i) = q_dot(i); %Adiciona a velocidade angular de flexão do elemento
    
    [~,~,~] = aeroposition_3d(n,i,L0,c,a); %Essa função calcula a posição
    % do centro aerodinâmico do elemento contando flexão e torção, agora
    % não estou usando ela, porque achei melhor separar o problema em
    % flexão e torção calculando eles separadamente na linha X.
    [~,Jpos,J_z] = position_elastic(n,i,L0); %Essa função calcula a posição
    %do centro elastico do elemento somente em função da flexão, isso
    %significa que é como o offset do centro elástico fosse 0 e não existe
    %torção atuando. Essa função devolve o jacobiano das coordenadas x,y,z
    %em função de q (Jpos) e J_z é o jacobiano da coordenada z em função de
    %q, é a velocidade vertical que usa no cálculo da sustentação.
    z_dot = J_z*qq_dot; %Velocidade vertical
    Bdot_T(n+i) = 1; % incremento no vetor de torção
    nf = [0;0;1];%[0;-sin(qq_T);cos(qq_T)]; %Direção da sustentação (vertical
    %para cima)
    
    % CALCULO DA AERODINAMICA LOCAL:
    % precisa do z_dot, qq_T, qq_Tdot local:
    L = lift_qs(U,rho,c,a,L0,qq_T,qq_Tdot,z_dot); %Cálculo da sustentação
    %no modelo quase-estacionário, eu já calibrei essa função com os
    %resultados da tese do Marko, então acho que esse não seja o problema
    M = M_qs(U,rho,L0,c,qq_Tdot);%Cálculo do segundo termo de momento
    %no modelo quase estacionário, já calibrei com o Marko também;
    
    % CALCULO DOS ESFORCOS GENERALIZADOS
    B1 = transpose(Jpos)*nf*L; % Termo do vetor de forças generalizadas
    %associado a flexão. Esse é o termo que acho que está gerando
    %problema, porque os outros dois a principio estão validados com a
    %dissertação do Marko. Ele calcula o jacobiano da posição global em
    %relação a q, multiplica pelo vetor de direção da força e multiplica
    %pela sustentação, acho que nesse jacobiano, deve estar entrando alguma
    %coisa a mais.
    B2 = Bdot_T*(0.25*c-a)*L; %Termo de torção associado a distância da
    %do centro aerodinamico para o centro elástico.
    B3 = - Bdot_T*M; % Segundo termo de momento no modelo quase estacionário
    B = B+ B1 + B2 + B3; %Vetor de saída, vetor de forças generalizadas
    L_vec(i) = L; %Vetor de sustentação para o plot
end

end