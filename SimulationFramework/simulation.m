function [X,Xdot] = simulation(beam_data,model,DoF,gravity,tspan,M,I,K,C,X0,varargin)

addpath background\utils
addpath background

global beam
beam = beam_data;
n_DoF = length(DoF)*sum(cat(1,beam.n_element));
%% Full-Order Model Simulation
if any(strcmp(model,'FOM')) || any(strcmp(model,'BOTH'))
    if ~any(strcmp(varargin,'PlotProgress'))
        opts = odeset('OutputFcn',@(t,X,flag) ode_progress(t,X,flag));
    else
        opts = odeset('OutputFcn',@(t,X,flag) ode_progress(t,X,flag,'PlotProgress'));
    end
    [T,X] = ode15s(@(t,X) dynamics(t,X,M,I,K,C,DoF,gravity,'FOM'),tspan,X0,opts);

    Xdot = zeros(2*n_DoF,length(T));
    for i = 1:length(T)
        Xdot(:,i) = dynamics(T(i),X(i,:).',M,I,K,C,DoF,gravity,'FOM');
    end
    Xdot = Xdot.';
end

end