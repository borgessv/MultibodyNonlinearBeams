function elements = create_elements(n_element,L)

syms q [4*n_element 1]
elements = struct("ID",{},"L",{},"C_ve",{},"r_v",{});
L_element = L/n_element;
r_aft_v = [0; 0; 0];

for i_element = 1:n_element
C_ev = DCM(1,q(i_element+n_element))*DCM(2,q(i_element+2*n_element))*DCM(3,q(i_element+3*n_element));
C_ve = C_ev.';

r_ea_e(:,2*i_element-1) = C_ev*r_aft_v;
r_ea_e(:,2*i_element) = r_ea_e(:,2*i_element-1) + [L_element + q(i_element); 0; 0];

r_ea_v(:,2*i_element-1) = C_ve*r_ea_e(:,2*i_element-1);
r_ea_v(:,2*i_element) = C_ve*r_ea_e(:,2*i_element);

r_aft_v = r_ea_v(:,2*i_element);

end
end