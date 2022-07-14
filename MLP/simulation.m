function [X,Xdot,H] = simulation(beam,model,DoF,gravity,tspan,X0,varargin)

addpath 'C:\Users\Vitor\Documents\Arquivos\Mestrado ITA\Dissertação\Codes\MultibodyNonlinearBeams'
addpath 'C:\Users\Vitor\Documents\Arquivos\Mestrado ITA\Dissertação\Codes\MultibodyNonlinearBeams\SupportCodes'

% Mass, stiffness and damping matrices:
progressbar('initializing matrices...')
%tspan = t0:dt:t1;
n_DoF = length(DoF)*sum(vertcat(beam.n_element));
[M,I] = mass_matrix(DoF);
K = stiffness_matrix(DoF);
%K(1,1) = 0;
c = 0;
C = c.*K;
if any(strcmp(DoF,'Torsion'))
    C(end-sum(beam.n_element)+1:end,end-sum(beam.n_element)+1:end) = 0.05*C(end-sum(beam.n_element)+1:end,end-sum(beam.n_element)+1:end);
end
progressbar('concluded')


%% Full-Order Model Simulation

if any(strcmp(model,'FOM')) || any(strcmp(model,'BOTH'))
    % FOM equilibrium solution:
%     progressbar('solving equilibrium...')
%     x0 = zeros(n_DoF,1);
%     options = optimoptions('fsolve','Display','off');
%     [xeq,~,exitflag] = fsolve(@(X) equilibrium(X,K,DoF,gravity,'FOM'),x0,options);
%     if any(exitflag,1:4)
%         progressbar('concluded')
%     else
%         error('Equilibrium solution could not be found or is not reliable!')
%     end

    % Plotting the FOM equilibrium condition:
%     if any(strcmp(model,'FOM'))
%         plot_structure(xeq,DoF,element)
%     end

    % Full-order model dynamics simulation:
    %X0 = zeros(2*n_DoF,1);
    %X0(n_DoF+1:end) = repmat(-linspace(pi/18,pi/12,n_DoF/length(DoF)),1,1).';
    %opts = odeset('OutputFcn',@(t,X,flag) ode_progress(t,X,flag,'PlotProgress'));
    [T,X] = ode15s(@(t,X) dynamics(t,X,M,I,K,C,DoF,gravity,'FOM'),tspan,X0);
    q_FOM = X(:,n_DoF+1:end).';
    p_FOM = X(:,1:n_DoF).';

    H = zeros(1,length(T));
    %     M_mat = kron(eye(length(DoF)),M(1:n_DoF/length(DoF),1:n_DoF/length(DoF)));
    %     if any(strcmp(DoF,'Torsion'))
    %         M_mat(n_DoF-n_DoF/length(DoF)+1:end,n_DoF-n_DoF/length(DoF)+1:end) = I(n_DoF-n_DoF/length(DoF)+1:end,n_DoF-n_DoF/length(DoF)+1:end);
    %     else
    %     end
    for i = 1:length(T)
        J = complexstep(@(q) element_positionCM(q,DoF),q_FOM(:,i));
        %         [~,J] = dlfeval(@(q) element_positionCM(q,DoF,'gradient'),dlarray(q_FOM(:,i)));
        %         J = extractdata(J);
        [~,~,~,Ug] = generalized_force(T(i),q_FOM(:,i),J,DoF,gravity);
        Xdot(:,i) = dynamics(T(i),X(i,:).',M,I,K,C,DoF,gravity,'FOM');
        %         qdot =  Xdot(n_DoF+1:end);
        H(i) = 0.5*(p_FOM(:,i).'*((J.'*M*J + I)\p_FOM(:,i)) + q_FOM(:,i).'*(K*q_FOM(:,i))) + Ug;
    end
end
end