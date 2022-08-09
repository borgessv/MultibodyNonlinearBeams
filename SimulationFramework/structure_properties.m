function [M,I,K,C] = structure_properties(beam_data,DoF,disp_progress)

global beam
beam = create_beam(beam_data,disp_progress);
n_element = sum(cat(1,beam.n_element));

[M,I] = mass_matrix(DoF,disp_progress);

K = stiffness_matrix(DoF,disp_progress);

c = 0.5;
C = damping_matrix(DoF,K,c,disp_progress);
if any(strcmp(DoF,'Torsion'))
    C(end-n_element:end) = 100*C(end-n_element:end);
end

filename = fullfile(sprintf('..\\SimulationFramework\\background\\beam_data.mat'));
save(filename,'beam')
end