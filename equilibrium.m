function Xeq = equilibrium(X,K,DoF)
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

Jq = lineariza_complexstep(@(q) element_position(q,DoF),X); % Numerical calculation of the Jacobian of the generalized coordinates using the complex step approach
Qe = generalized_force(0,Jq,DoF,'GravityOn'); % Calculation of generalized forces acting on the system

Xeq = K*X - Qe;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%