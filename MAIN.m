clear all % 'clear all' must be used since persistent variables might have been used in a previous run
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 'FOM'; % Options: 'FOM' for full-order model, 'ROM' for reduced order model by truncation method or 'BOTH' for both models analysis in a single run
DoF = {'OutBend','Torsion','InBend'}; % Options: 'InBend': in-plane bending; 'OutBend': out-of-plane bending; 'Axial': axial deformation; 'Torsion': torsion
gravity = 'GravityOn'; % Options: 'GravityOn' to consider gravitational force or 'GravityOff' to disconsider it. 
n_modes = 10; % Number of modes used to produce a ROM using truncation method
t0 = 0; % Initial time of the simulation [s]
dt = 0.1; % Time step of simulation [s]
t1 = 40; % End time of simulation [s]
animate = 'Y'; % Options: 'Y' or 'N'
elem_type = '3D'; % Options: '1D' for unidimensional elements or '3D' for 3-dimensional elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------<strong>ANALYSIS IN PROGRESS</strong>----------------------------------\n')
addpath 'Support Codes'
global beam
load beams_test.mat beam % The .mat file must be created previously in 'create_beam.m'

if any(strcmp(model,'FOM'))
    fprintf('<strong>Model</strong>: Full-Order \n<strong>Degrees of Freedom</strong>: ')
    fprintf('%s, ',DoF{1:end-1})
    fprintf('%s\n',DoF{end})
    fprintf('<strong>Number of elements</strong>: %d\n',sum(beam.n_element))
elseif any(strcmp(model,'ROM'))
    fprintf('<strong>Model</strong>: Reduced-Order \n<strong>Degrees of Freedom</strong>: ')
    fprintf('%s, ',DoF{1:end-1})
    fprintf('%s\n',DoF{end})
    fprintf('<strong>Number of elements</strong>: %d\n',sum(beam.n_element))
elseif any(strcmp(model,'BOTH'))
    fprintf('<strong>Model</strong>: Full-Order and Reduced-Order \n<strong>Degrees of Freedom</strong>: ')
    fprintf('%s, ',DoF{1:end-1})
    fprintf('%s\n',DoF{end})
    fprintf('<strong>Number of elements</strong>: %d\n',sum(beam.n_element))
end

% Mass, stiffness and damping matrices:
progressbar('initializing matrices...')
tspan = t0:dt:t1;
n_DoF = length(DoF)*sum(vertcat(beam.n_element));
M = mass_matrix(DoF);
K = stiffness_matrix(DoF);
c = 0.1;
C = c.*K;
C(end-sum(beam.n_element)+1:end,end-sum(beam.n_element)+1:end) = 0.005*C(end-sum(beam.n_element)+1:end,end-sum(beam.n_element)+1:end);
progressbar('concluded')


%% Full-Order Model Simulation

if any(strcmp(model,'FOM')) || any(strcmp(model,'BOTH'))
    % FOM equilibrium solution:
    progressbar('solving equilibrium...')
    x0 = zeros(n_DoF,1);
    options = optimoptions('fsolve','Display','off');
    [xeq,~,exitflag] = fsolve(@(X) equilibrium(X,K,DoF,gravity,'FOM'),x0,options);
    if any(exitflag,1:4)
        progressbar('concluded')
    else
        error('Equilibrium solution could not be found or is not reliable!')
    end

    % Plotting the FOM equilibrium condition:
    if any(strcmp(model,'FOM'))
        plot_structure(xeq,DoF,elem_type)
    end

    % Full-order model dynamics simulation:
    X0 = zeros(2*n_DoF,1);
    X0(n_DoF+1:end) = xeq;
    opts = odeset('OutputFcn',@(t,X,flag) ode_progress(t,X,flag,'PlotProgress'));
    [T,X] = ode15s(@(t,X) dynamics(t,X,M,K,C,DoF,gravity,'FOM'),tspan,X0,opts);
    q_FOM = X(:,n_DoF+1:end);

%     for i = 1:length(T)
%         Xdot(i,:) = dynamics(T(i),X(i,:).',M,K,C,DoF,gravity,'FOM');
%     end
% 
%     [x,y] = meshgrid(linspace(-5,5,401),linspace(-25,25,401));


    % Create animation of the FOM dynamics solution:
    if any(strcmp(model,'FOM')) && any(strcmp(animate,'Y'))
        animation(T,q_FOM,DoF,elem_type,'FOM_sim_test_followerforce','gif','avi')
    end
end


%% Reduced Order Model Simulation - Truncation Method

if any(strcmp(model,'ROM')) || any(strcmp(model,'BOTH'))
    % Initialization of the reduced order model:
    progressbar('creating reduced order model...')
    if any(strcmp(model,'ROM'))
        x0 = zeros(n_DoF,1);
        options = optimoptions('fsolve','Display','off');
        [xeq,~,exitflag] = fsolve(@(X) equilibrium(X,K,DoF,gravity,'FOM'),x0,options);
        if ~any(exitflag,1:4)
            error('Equilibrium solution could not be found or is not reliable!')
        end
    end
    Jq_eq = lineariza_complexstep(@(q) element_positionCM(q,DoF),xeq);
    I_i = 0.1;
    T = blkdiag(zeros(n_DoF-sum(vertcat(beam.n_element))),eye(sum(vertcat(beam.n_element))));
    I = I_i*eye(n_DoF);
    [eigvec,eigval] = eig(K,(Jq_eq.'*M*Jq_eq + T.'*I*T));
    lambda = diag(eigval);
%     lambda = reshape(lambda,[sum(vertcat(beam.n_element)),length(DoF)]);
%     [eigval,i_order] = sort(lambda,1);
%     i_order = i_order + (0:sum(vertcat(beam.n_element)):(length(DoF)-1)*sum(vertcat(beam.n_element)));
%     i_order = i_order.';
    [eigval,i_order] = sort(lambda);
    phi = eigvec(:,i_order);
    phi_r = phi(:,1:n_modes);
    progressbar('concluded')

    % ROM equilibrium solution:
    progressbar('solving ROM equilibrium...')
    eta0 = zeros(n_modes,1);
    options = optimoptions('fsolve','Display','off');
    [eta_eq,~,exitflag] = fsolve(@(X) equilibrium(X,K,DoF,gravity,'ROM',phi_r),eta0,options);
    if any(exitflag,1:4)
        progressbar('concluded')
        xeq_r = phi_r*eta_eq;
    else
        error('Equilibrium solution could not be found or is not reliable!')
    end

    % Plotting the ROM equilibrium condition:
    if any(strcmp(model,'ROM'))
        plot_structure(xeq_r,DoF,elem_type)
        drawnow
    elseif any(strcmp(model,'BOTH'))
        plot_structure(xeq,DoF,elem_type,xeq_r)
        drawnow
    end

    % Reduced-order model dynamics simulation:
    X0_r = zeros(2*n_modes,1);
    X0_r(n_modes+1:end) = eta_eq;
    opts = odeset('OutputFcn',@(t,eta,flag) ode_progress(t,eta,flag,'PlotProgress'));
    [T,X_r] = ode15s(@(t,eta) dynamics(t,eta,M,K,C,DoF,gravity,'ROM',phi_r),tspan,X0_r,opts);
    q_ROM = (phi_r*X_r(:,n_modes+1:end).').';

    % Create animation of the ROM dynamics solution:
    if any(strcmp(model,'ROM')) && any(strcmp(animate,'Y'))
        animation(T,q_ROM,DoF,elem_type,'ROMsim_test_followerforce','gif','avi')
    end
end

if any(strcmp(model,'BOTH')) && any(strcmp(animate,'Y'))
    animation(T,q_FOM,DoF,elem_type,'test_50_S1223','gif','avi',q_ROM)
end

fprintf('-----------------------------------<strong>END OF ANALYSIS</strong>------------------------------------\n')