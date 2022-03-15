function estrutura = goland
    L = 20*0.3048;        % [m]      
    b = 3*0.3048;         % [ft]
    m = 0.746*47.880259;     % [kg/m]
    GJ = 2.39e6*0.413253311;   % [lbf.ft²]
    EI = 23.65e6*0.413253311;  % [lbf.ft²]
    ra = 0.5; 
    xa = 0.1;
    I = m*(b*ra)^2; % kg m^2 (inercia secao torcao)
    a = -0.2;
    estrutura = struct('L',L,'I',I,'xa',xa,'b',b,'m',m,'GJ',GJ,'EI',EI, 'a', a);
end