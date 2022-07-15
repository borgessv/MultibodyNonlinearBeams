clear all; close all; clc; % 'clear all' must be used since persistent variables might have been used in a previous run
fprintf('--------------------------------<strong>ANALYSIS IN PROGRESS</strong>----------------------------------\n')
addpath background\utils
addpath background

global beam
beam = create_beam('beam_data_test.xlsx');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL DETAILS:
model = 'FOM'; % Options: 'FOM' for full-order model, 'ROM' for reduced order model by truncation method or 'BOTH' for both models analysis in a single run
DoF = {'OutBend','Torsion'}; % Options: 'InBend': in-plane bending; 'OutBend': out-of-plane bending; 'Axial': axial deformation; 'Torsion': torsion
gravity = 'GravityOn'; % Options: 'GravityOn' to consider gravitational force or 'GravityOff' to disconsider it.
tspan = 0:0.1:30; % Period of simulation [s]

% INITIAL CONDITIONS - 'FOM' OR 'BOTH' MODEL:
n_DoF = length(DoF)*sum(cat(1,beam.n_element)); 
p0 = zeros(n_DoF,1); % Generalized momenta initial condition
q0 = 'Xeq'; % Generalized coordinates initial condition

% INITIAL CONDITIONS - 'ROM' OR 'BOTH' MODEL:
n_modes = 1; % Number of modes used to produce a ROM using truncation method
eta_p0 = zeros(n_modes,1); % Generalized modal momenta initial conditions
eta_q0 = 'Xeq'; % Generalized modal coordinates initial conditions

% ANIMATION SETTINGS:
animate = 'Y'; % Options: 'Y' or 'N'
element = '3D'; % Options: '1D' for unidimensional elements or '3D' for 3-dimensional elements
animation_file = 'sim_test'; % Filename for the simulation [str]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializing Matrices:
[M,I] = mass_matrix(DoF);
K = stiffness_matrix(DoF);
C = damping_matrix(DoF,K,0.1);

% Equilibrium Solution:
if any(strcmp(q0,'Xeq')) || any(strcmp(model,'ROM')) || any(strcmp(model,'BOTH'))
    x0 = zeros(n_DoF,1);
    Xeq = solve_equilibrium(K,DoF,gravity,model,x0,'PlotEq',element);
    q0 = Xeq;
end

% FOM Simulation:
if any(strcmp(model,'FOM')) || any(strcmp(model,'BOTH'))
    X0 = [p0;q0];
    [X,Xdot] = simulation(beam,model,DoF,gravity,tspan,M,I,K,C,X0);
    p_FOM = X(:,1:n_DoF);
    q_FOM = X(:,n_DoF+1:end);

    % Hamiltonian:
    H = zeros(length(tspan));
    for i = 1:length(tspan)
        J = complexstep(@(q) element_positionCM(q,DoF),q_FOM(i,:));
        [~,~,~,Ug] = generalized_force(tspan(i),q_FOM(i,:),J,DoF,gravity);
        H(i) = 0.5*(p_FOM(i,:)*((J.'*M*J + I)\p_FOM(i,:).') + q_FOM(i,:)*(K*q_FOM(i,:).')) + Ug;
    end

    % Simulation Results Animation:
    if any(strcmp(model,'FOM')) && any(strcmp(animate,'Y'))
        animation(tspan,q_FOM,H,DoF,element,[animation_file,'_FOM'],'gif','avi')
    end
end


%% Reduced Order Model Simulation - Truncation Method

if any(strcmp(model,'ROM')) || any(strcmp(model,'BOTH'))
    % Initialization of the reduced order model:
    phi_r = create_rom(n_modes,DoF,M,I,K,Xeq) ;

    % Equilibrium Solution:
    if any(strcmp(eta_q0,'Xeq'))
        eta0 = zeros(n_modes,1);
        eta_eq = solve_equilibrium(K,DoF,gravity,model,eta0,phi_r,Xeq,'PlotEq',element);
        eta_q0 = eta_eq;
    end

    % Reduced-order model dynamics simulation:
    X0_ROM = [eta_p0,eta_q0];
    opts = odeset('OutputFcn',@(t,eta,flag) ode_progress(t,eta,flag,'PlotProgress'));
    [T,X_r] = ode15s(@(t,eta) dynamics(t,eta,M,I,K,C,DoF,gravity,'ROM',phi_r),tspan,X0_ROM,opts);
    q_ROM = phi_r*X_r(:,n_modes+1:end).';
    p_ROM = phi_r*X_r(:,1:n_modes).';

    H_ROM = zeros(1,length(T));
    M_mat = kron(eye(length(DoF)),M(1:n_DoF/length(DoF),1:n_DoF/length(DoF)));
    if any(strcmp(DoF,'Torsion'))
        M_mat(n_DoF-n_DoF/length(DoF)+1:end,n_DoF-n_DoF/length(DoF)+1:end) = I(n_DoF-n_DoF/length(DoF)+1:end,n_DoF-n_DoF/length(DoF)+1:end);
    else
    end
    for i = 1:length(T)
        H_ROM(i) = 0.5*(p_ROM(:,i).'*(M_mat\p_ROM(:,i)) + q_ROM(:,i).'*(K*q_ROM(:,i)));
    end

    % Create animation of the ROM dynamics solution:
    if any(strcmp(model,'ROM')) && any(strcmp(animate,'Y'))
        animation(T,q_ROM.',H_ROM,DoF,element,[animation_file,'_ROM'],'gif','avi')
    end
end

if any(strcmp(model,'BOTH')) && any(strcmp(animate,'Y'))
    animation(T,q_FOM.',H,DoF,element,[animation_file,'_BOTH'],'gif','avi',q_ROM.',H_ROM)
end

fprintf('-----------------------------------<strong>END OF ANALYSIS</strong>------------------------------------\n')