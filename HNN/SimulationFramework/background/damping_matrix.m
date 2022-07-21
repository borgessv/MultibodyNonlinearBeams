function C = damping_matrix(DoF,varargin)
progressbar('creating damping matrix...')

global beam

K = varargin{1};
c = varargin{2};

C = c.*K;
if any(strcmp(DoF,'Torsion'))
    C(end-sum(beam.n_element)+1:end,end-sum(beam.n_element)+1:end) = 0.05*C(end-sum(beam.n_element)+1:end,end-sum(beam.n_element)+1:end);
end
progressbar('done')

end