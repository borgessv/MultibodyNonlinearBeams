function estrutura = highlyflex

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

    L = 16;        % [m]      
    b = 1/2;         % [m]
    m = 0.75;     % [kg/m]
    GJ = 1e4;   % [N.m²]
    EI = 2e4;  % [N.m²]
    I = 0.1; 
    xa = 0;
    a = 0;
    estrutura = struct('L',L,'I',I,'xa',xa,'b',b,'m',m,'GJ',GJ,'EI',EI, 'a', a);
end