function [M,I,K,C] = structure_properties(beam_data,DoF)
addpath background
global beam

beam = create_beam(beam_data);
[M,I] = mass_matrix(DoF);
K = stiffness_matrix(DoF);
c = 0.1;
C = damping_matrix(DoF,K,c);

filepath = 'C:\Users\Vitor\Documents\Arquivos\MastersDegree\MSc_Thesis\Codes\MultibodyNonlinearBeams\SimulationFramework\background';
filename = fullfile(filepath, sprintf('beam_data.mat'));
save(filename,'beam')
end