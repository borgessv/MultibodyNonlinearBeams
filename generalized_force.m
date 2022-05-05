function Qe = generalized_force(t,Jq,DoF,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION OVERVIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The function generalized_force(t,Jq,varargin) returns the generalized 
%   forces acting upon the system. 
%   IMPORTANT: If external forces besides the gravitational one are present 
%   their vectors must be specified - see lines 48-80.
%   INPUTS:
%       t: time;
%       Jq: Jacobian matrix evaluated previously using e.g. complex step;
%       varargin*: specifies wheter or not to consider gravity force.
%                  'GravityOn': cosiders gravity force;
%                  'GravityOff': does not consider gravity force;
%       *if not declared, the 'GravityOn' option is assumed.  
%   OUTPUT:
%       - Qe: Column vector of generalized forces acting upon the system,
%       such that each element is the generalized force related to the
%       correspondent degree of freedom equation.
%
%   Author: Vitor Borges Santos - borgessv93@gmail.com
%   Date: 15/04/22
%
%%%%%%%%%%%%%%%%%%%%%%%% BEGINNING OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
global beam
% Generalized forces due to gravity:
if any(strcmp(varargin,'GravityOn')) || isempty(varargin)
    g = 9.80665; % [m/s²] - gravitational acceleration assuming a constant gravitational field at sea level and geodetic latitude of 45°
    for i_beam = 1:length(beam)
        n = beam(i_beam).n_element;
        Fg_vec = zeros(3,n);
        Qg_element = zeros(n,length(DoF)*n);
        for i_element = 1:n
            Fg_vec(:,i_element) = [0; 0; -beam(i_beam).element(i_element).m*g];
            for i_q = 1:length(DoF)*n
                Qg_element(i_element,i_q) = dot(Fg_vec(:,i_element),[Jq(i_element,i_q); Jq(n+i_element,i_q); Jq(2*n+i_element,i_q)]);
            end
        end
        beam(i_beam).Qg = sum(Qg_element,1);
    end
    Qg = horzcat(beam.Qg).';
elseif any(strcmp(varargin,'GravityOff'))
    for i_beam = 1:length(beam)
        n = beam(i_beam).n_element;
        beam(i_beam).Qg = zeros(1,length(DoF)*n);
    end
    Qg = horzcat(beam.Qg).';
end

% Generalized forces due to other external forces:
n_beam = length(beam);
for i_beam = 1:n_beam
    n = beam(i_beam).n_element;
    F_vec = zeros(3,n);
    Q_element = zeros(n,length(DoF)*n);
    %%%%%%%%%%%%% Repeat this block to consider different forces acting on different beams %%%%%%%%%%%%%
    if any(i_beam == 1)
        for i_element = 1:n
            %%%%%%%% Repeat this block to consider different forces acting on different elements %%%%%%%
            if any(i_element == n)
                if t >= 2 && t <= 4
                    F_vec(:,i_element) = 500*[0; -10*sin(pi/2*t); 2*sin(pi/2*t)]; % External force acting on each element of the beam
                elseif t>0 && t < 2
                    F_vec(:,i_element) = [0; -5000*sin(pi/2*t); 0];
                else
                    F_vec(:,i_element) = [-0; -0; -0];
                end
                for i_q = 1:length(DoF)*n
                    Q_element(i_element,i_q) = dot(F_vec(:,i_element),[Jq(i_element,i_q); Jq(n+i_element,i_q); Jq(2*n+i_element,i_q)]);
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of element block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                Q_element(i_element,1:length(DoF)*n) = zeros(1,length(DoF)*n);
            end
        end
        beam(i_beam).Q = sum(Q_element,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of beam block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        beam(i_beam).Q = zeros(1,length(DoF)*n);
    end
end
Q = horzcat(beam.Q).';

% Total generalized forces acting on the system:
Qe = Qg + Q;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%