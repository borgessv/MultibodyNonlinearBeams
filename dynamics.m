function Xdot = dynamics(t,X,M,K,C,DoF)
%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION OVERVIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The function dynamics(t,X,M,K,C) returns the equations of motion of the
%   system using Hamiltonian mechanics. 
%   INPUTS: 
%       t: time;
%       X: generalized coordinates vector;
%       M: mass matrix of the system (obtained previously by the function
%       mass_matrix);
%       K: stiffness matrix of the system (obtained previously by the 
%       function stiffness_matrix);
%       C: damping matrix (must be calculated previously or assumed e.g. as
%       c*K, where c is the damping coefficient.
%   OUTPUT:
%       Xdot: vector containing the time derivatives of the generalized 
%       momentum and coordinates, i.e. the Hamiltonian equations of motion.  
%
%   Author: Vitor Borges Santos - borgessv93@gmail.com
%   Date: 15/04/22
%
%%%%%%%%%%%%%%%%%%%%%%%% BEGINNING OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
global beam
n_DoF = length(DoF)*sum(vertcat(beam.n_element));
Jq = lineariza_complexstep(@(q) element_position(q,DoF),X(n_DoF+1:end)); % Numerical calculation of the Jacobian of the generalized coordinates using the complex step approach
Q = generalized_force(t,Jq,DoF,'GravityOn');

qdot = (Jq.'*M*Jq)\X(1:n_DoF);
pdot = -K*X(n_DoF+1:end) + Q - C*qdot;

Xdot = [pdot;qdot];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%