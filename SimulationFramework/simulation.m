function [X,Xdot] = simulation(model,DoF,gravity,tspan,M,I,K,C,X0,varargin)

addpath background
addpath background\CrossSectionData
addpath background\utils

global beam
load beam_data.mat beam
n_DoF = length(DoF)*sum(cat(1,beam.n_element));
if ~any(strcmp(varargin,'PlotProgress'))
    opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'OutputFcn',@(t,X,flag) ode_progress(t,X,flag));
else
    opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'OutputFcn',@(t,X,flag) ode_progress(t,X,flag,'PlotProgress'));
end

%% Full-Order Model Simulation
if any(strcmp(model,'FOM')) || (any(strcmp(model,'BOTH')) && ~any(cellfun(@isnumeric,varargin)))

    [T,X] = ode15s(@(t,X) dynamics(t,X,M,I,K,C,DoF,gravity,'FOM'),tspan,X0,opts);

    Xdot = zeros(2*n_DoF,length(T));
    for i = 1:length(T)
        Xdot(:,i) = dynamics(T(i),X(i,:).',M,I,K,C,DoF,gravity,'FOM');
    end
    Xdot = Xdot.';

%% Reduced-Order Model Simulation:
elseif any(strcmp(model,'ROM')) || (any(strcmp(model,'BOTH')) && any(cellfun(@isnumeric,varargin)))

    phi_r = varargin{1};
    [T,X] = ode15s(@(t,eta) dynamics(t,eta,M,I,K,C,DoF,gravity,'ROM',phi_r),tspan,X0,opts);

    Xdot = zeros(size(X,2),length(T));
    for i = 1:length(T)
        Xdot(:,i) = dynamics(T(i),X(i,:).',M,I,K,C,DoF,gravity,'ROM',phi_r);
    end
    Xdot = Xdot.';
end

end