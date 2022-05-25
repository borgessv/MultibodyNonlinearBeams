function xeq = equilibrium(X,K,DoF,gravity,model,phi_r)
%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION OVERVIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The function equilibrium(X,K) returns the equilibrium solution of the
%   at a given condition.
%   INPUTS:
%       X: generalized coordinates vector;
%       K: stiffness matrix of the system (obtained by the function
%       stiffness_matrix).
%   OUTPUT:
%       Xeq: vector containing the values of the generalized coordinates
%       obtained in the equilibrium solution.
%
%   Author: Vitor Borges Santos - borgessv93@gmail.com
%   Date: 15/04/22
%
%%%%%%%%%%%%%%%%%%%%%%%% BEGINNING OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(model,'FOM'))
    Jq = lineariza_complexstep(@(q) element_position(q,DoF),X); % Numerical calculation of the Jacobian of the generalized coordinates using the complex step approach
    Qe = generalized_force(0,X,Jq,DoF,gravity); % Calculation of generalized forces acting on the system

    xeq = K*X - Qe;

elseif any(strcmp(model,'ROM'))
    Jq = lineariza_complexstep(@(q) element_position(q,DoF),phi_r*X); % Numerical calculation of the Jacobian of the generalized coordinates using the complex step approach
    Qe = generalized_force(0,phi_r*X,Jq,DoF,gravity); % Calculation of generalized forces acting on the system

    xeq = (phi_r.'*K*phi_r)*X - phi_r.'*Qe;

else
    error('The type of the model must be specified! Choose "FOM" if full-order model or "ROM" if reduced order model.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%