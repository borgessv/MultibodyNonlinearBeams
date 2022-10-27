function phi_r = create_rom(n_modes,DoF,Xeq)
progressbar('creating reduced order model...')
load beam_data.mat M I K
Jq_eq = complexstep(@(q) element_positionCM(q,DoF),Xeq);
[eigvec,eigval] = eig(K,(Jq_eq.'*M*Jq_eq + I));
lambda = diag(eigval);
%     lambda = reshape(lambda,[sum(vertcat(beam.n_element)),length(DoF)]);
%     [eigval,i_order] = sort(lambda,1);
%     i_order = i_order + (0:sum(vertcat(beam.n_element)):(length(DoF)-1)*sum(vertcat(beam.n_element)));
%     i_order = i_order.';
[~,i_order] = sort(lambda);
phi = eigvec(:,i_order);
phi_r = phi(:,1:n_modes);
progressbar('done')
end