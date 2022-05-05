function M = mass_matrix(DoF)
%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION OVERVIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The function mass_matrix returns the mass matrix of the system. 
%   INPUTS: N/A
%   OUTPUT:
%       M: diagonal matrix containing the masses of each element of each 
%       beam.
%
%   Author: Vitor Borges Santos - borgessv93@gmail.com
%   Date: 15/04/22
%
%%%%%%%%%%%%%%%%%%%%%%%% BEGINNING OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
global beam
n_beam = length(beam);
for i_beam = 1:n_beam
    n = beam(i_beam).n_element;
    M_element = zeros(1,n);
    for i_element = 1:n
        M_element(i_element) = beam(i_beam).element(i_element).m;
    end
    beam(i_beam).M = diag(repmat(M_element,1,3));
end
M = blkdiag(beam.M);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%