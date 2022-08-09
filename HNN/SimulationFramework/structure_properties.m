function [M,I,K,C] = structure_properties(beam_data,DoF)
addpath background
global beam

beam = create_beam(beam_data);
[M,I] = mass_matrix(DoF);
K = stiffness_matrix(DoF);
c = 0.00005;
C = damping_matrix(DoF,K,c);

%filepath = fileparts(cd);
filename = fullfile(sprintf('..\\SimulationFramework\\background\\beam_data.mat'));
save(filename,'beam')
end