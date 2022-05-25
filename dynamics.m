function Xdot = dynamics(t,X,M,K,C,DoF,gravity,model,phi_r)
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
if any(strcmp(model,'FOM'))
    global beam
    n_DoF = length(DoF)*sum(vertcat(beam.n_element));
    Jq = lineariza_complexstep(@(q) element_position(q,DoF),X(n_DoF+1:end)); % Numerical calculation of the Jacobian of the generalized coordinates using the complex step approach
    Q = generalized_force(t,X(n_DoF+1:end),Jq,DoF,gravity);

    qdot = (Jq.'*M*Jq)\X(1:n_DoF);
    pdot = -K*X(n_DoF+1:end) + Q - C*qdot;

    Xdot = [pdot;qdot];

elseif any(strcmp(model,'ROM'))
    r = size(phi_r,2);
    Jq = lineariza_complexstep(@(q) element_position(q,DoF),phi_r*X(r+1:end)); % Numerical calculation of the Jacobian of the generalized coordinates using the complex step approach
    Q = generalized_force(t,phi_r*X(r+1:end),Jq,DoF,gravity);

    etadot = (phi_r.'*Jq.'*M*Jq*phi_r)\X(1:r);
    prdot = -(phi_r.'*K*phi_r)*X(r+1:end) + phi_r.'*Q - (phi_r.'*C*phi_r)*etadot;

    Xdot = [prdot;etadot];
else
    error('The type of the model must be specified! Choose "FOM" if full-order model or "ROM" if reduced order model')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%