clear all % 'clear all' must be used since persistent variables might have been used in a previous run
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 'FOM'; % Options: 'FOM' for full-order model, 'ROM' for reduced order model by truncation method or 'BOTH' for both models analysis in a single run
DoF = {'OutBend'}; % Options: 'InBend': in-plane bending; 'OutBend': out-of-plane bending; 'Axial': axial deformation; 'Torsion': torsion
gravity = 'GravityOn'; % Options: 'GravityOn' to consider gravitational force or 'GravityOff' to disconsider it. 
n_modes = 30; % Number of modes used to produce a ROM using truncation method
t0 = 0; % Initial time of the simulation [s]
dt = 0.005; % Time step of simulation [s]
t1 = 120; % End time of simulation [s]
animate = 'Y'; % Options: 'Y' or 'N'
element = '3D'; % Options: '1D' for unidimensional elements or '3D' for 3-dimensional elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------<strong>ANALYSIS IN PROGRESS</strong>----------------------------------\n')
addpath 'SupportCodes'
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
[M,I] = mass_matrix(DoF);
K = stiffness_matrix(DoF);
%K(1,1)=0;
c = 0.1;
C = c.*K;
%K(1,1) = 0;
if any(strcmp(DoF,'Torsion'))
    C(end-sum(beam.n_element)+1:end,end-sum(beam.n_element)+1:end) = 0.005*C(end-sum(beam.n_element)+1:end,end-sum(beam.n_element)+1:end);
end
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
        plot_structure(xeq,DoF,element)
    end

    % Full-order model dynamics simulation:
    [p0,q0] = meshgrid(linspace(0,0,1),linspace(0,0,1));
    points = [p0(:),q0(:)];
    for k = 1:length(points)
    X0 = [0;xeq];%zeros(2*n_DoF,1);
%     X0 = [points(k,1);points(k,2)];
    %X0 = [pi*(2*randn(n_DoF,1)-1);pi/2*(2*randn(n_DoF,1)-1)];
    %X0(n_DoF+1:end) = repmat(-linspace(pi/18,pi/12,n_DoF/length(DoF)),1,1).';
    opts = odeset('OutputFcn',@(t,X,flag) ode_progress(t,X,flag,'PlotProgress'));
    [T,X] = ode15s(@(t,X) dynamics(t,X,M,I,K,C,DoF,gravity,'FOM'),tspan,X0);

    for i = 1:length(T)
        %J = complexstep(@(q) element_positionCM(q,DoF),q_FOM(:,i));
        %[~,~,Qext,Ug] = generalized_force(T(i),q_FOM(:,i),J,DoF,gravity);
         Xdot(:,i) = dynamics(T(i),X(i,:).',M,I,K,C,DoF,gravity,'FOM');         
         %         qdot =  Xdot(n_DoF+1:end);
        %H(i) = 0.5*(p_FOM(:,i).'*((J.'*M*J + I)\p_FOM(:,i)) + q_FOM(:,i).'*(K*q_FOM(:,i))) + Ug;
    end
        pdot_FOM(k,:) = Xdot(1:n_DoF,:).';
        qdot_FOM(k,:) = Xdot(n_DoF+1:end,:).';
    end
    %figure; plot(T,H)
    % Create animation of the FOM dynamics solution:
    if any(strcmp(model,'FOM')) && any(strcmp(animate,'Y'))
        animation(T,q_FOM.',H,DoF,element,'FOM_sim_test_followerforce','gif','avi')
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
    Jq_eq = complexstep(@(q) element_positionCM(q,DoF),xeq);
    [eigvec,eigval] = eig(K,(Jq_eq.'*M*Jq_eq + I));
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
        plot_structure(xeq_r,DoF,element)
        drawnow
    elseif any(strcmp(model,'BOTH'))
        plot_structure(xeq,DoF,element,xeq_r)
        drawnow
    end

    % Reduced-order model dynamics simulation:
    X0_r = zeros(2*n_modes,1);
    %X0_r(n_modes+1:end) = eta_eq;
    opts = odeset('OutputFcn',@(t,eta,flag) ode_progress(t,eta,flag,'PlotProgress'));
    [T,X_r] = ode15s(@(t,eta) dynamics(t,eta,M,I,K,C,DoF,gravity,'ROM',phi_r),tspan,X0_r,opts);
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
        animation(T,q_ROM.',H_ROM,DoF,element,'ROMsim_test_followerforce','gif','avi')
    end
end

if any(strcmp(model,'BOTH')) && any(strcmp(animate,'Y'))
    animation(T,q_FOM.',H,DoF,element,'test_50_S1223','gif','avi',q_ROM.',H_ROM)
end

fprintf('-----------------------------------<strong>END OF ANALYSIS</strong>------------------------------------\n')