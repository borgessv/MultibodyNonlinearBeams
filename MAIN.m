clear all % Must use 'clear all' to clear persistent variables that might have been used in a previous run
close all
clc
addpath 'Support Codes'

fprintf('---------------------------SIMULATION IN PROGRESS----------------------------\n')
global beam
load beams_test.mat beam

DoF = {'OutBend','InBend'}; % 'InBend': in-plane bending; 'OutBend': out-of-plane bending; 'Axial': axial deformation; 'Torsion': torsion
n_DoF = length(DoF)*sum(vertcat(beam.n_element));

% Mass, stiffness and damping matrices:
progressbar('initializng matrices...')
M = mass_matrix(DoF);
K = stiffness_matrix(DoF);
c = 0.1;
C = c.*K;
progressbar('concluded')

% Static solution - Equilibrium:
progressbar('solving equilibrium...')
x0 = zeros(n_DoF,1);
options = optimoptions('fsolve','Display','off');
[xeq,fval,exitflag] = fsolve(@(X) equilibrium(X,K,DoF),x0,options);
if exitflag == 1
    progressbar('concluded')
else
    error('Equilibrium solution could not be found or is not reliable!')
end

% Plotting the equilibrium:
r_eq = element_position(xeq,DoF);
r = [beam(1).r0(1,:);beam(1).r1];
rCM = beam(1).rCM;
figure
plot3(r(:,1),r(:,2),r(:,3),'-|k','markersize',2,'markerfacecolor','k')
hold on
plot3(rCM(:,1),rCM(:,2),rCM(:,3),'sr','markersize',3,'markerfacecolor','r')
title('Equilibrium Solution','interpreter','latex')
xlabel('x [m]','interpreter','latex')
ylabel('y [m]','interpreter','latex')
zlabel('z [m]','interpreter','latex')
grid on
axis equal
set(gcf,'color','w');
drawnow

% Dynamics solution:
X0 = zeros(2*n_DoF,1);
X0(n_DoF+1:end) = xeq;
t0 = 0;
dt = 0.1;
t1 = 25;
tspan = t0:dt:t1;
figure
opts = odeset('OutputFcn',@(t,X,flag) ode_progress(t,X,flag,'Plot'));
[T,X] = ode15s(@(t,X) dynamics(t,X,M,K,C,DoF),tspan,X0,opts);

% Create animation of the dynamics solution:
animation_dynamics(T,X,DoF,'sim_test','gif','avi')

fprintf('------------------------------END OF SIMULATION------------------------------\n')