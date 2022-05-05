function rCM = element_position(q,DoF)

global beam

for i_beam = 1:length(beam)
    n = beam(i_beam).n_element;
    %     if i_beam == 1
    %         q = X(1:3*n);
    %     else
    %         n_aux = sum(vertcat(beam(1:i_beam).n_element));
    %         q = X(3*(n_aux-n)+1:3*(n_aux-n)+3*n);
    %     end
    %     r0x = zeros(1,n);
    %     r0y = zeros(1,n);
    %     r0z = zeros(1,n);
    %     r1x = zeros(1,n);
    %     r1y = zeros(1,n);
    %     r1z = zeros(1,n);
    %     rCMx = zeros(1,n);
    %     rCMy = zeros(1,n);
    %     rCMz = zeros(1,n);
    beam(i_beam).r0 = zeros(n,3);
    beam(i_beam).r1 = zeros(n,3);
    beam(i_beam).rCM = zeros(n,3);
    r0_0_element = beam(i_beam).element(1).r0;
    for i_element = 1:n
        C_i0 = beam(i_beam).C_i0;
        C_di = eye(3);
        
        if any(strcmp(DoF,'InBend')) && any(strcmp(DoF,'OutBend'))
            C_di = DCM(3,sum(q(n+1:i_element+n)));
        elseif any(strcmp(DoF,'InBend')) && ~any(strcmp(DoF,'OutBend'))
            C_di = DCM(3,sum(q(1:i_element)));
        end
        
        if any(strcmp(DoF,'OutBend'))
            C_di = DCM(2,sum(q(1:i_element)))*C_di;
        end

        if any(strcmp(DoF,'Torsion')) && any(strcmp(DoF,'InBend')) && any(strcmp(DoF,'OutBend')) && any(strcmp(DoF,'Axial'))
            C_di = DCM(1,sum(q(3*n+1:i_element+3*n)))*C_di;
        elseif any(strcmp(DoF,'Torsion')) && any(strcmp(DoF,'InBend')) && any(strcmp(DoF,'OutBend')) && ~any(strcmp(DoF,'Axial')) || any(strcmp(DoF,'Torsion')) && any(strcmp(DoF,'InBend')) && ~any(strcmp(DoF,'OutBend')) && any(strcmp(DoF,'Axial')) || any(strcmp(DoF,'Torsion')) && ~any(strcmp(DoF,'InBend')) && any(strcmp(DoF,'OutBend')) && any(strcmp(DoF,'Axial'))
            C_di = DCM(1,sum(q(2*n+1:i_element+2*n)))*C_di;
        elseif any(strcmp(DoF,'Torsion')) && ~any(strcmp(DoF,'InBend')) && any(strcmp(DoF,'OutBend')) && ~any(strcmp(DoF,'Axial')) || any(strcmp(DoF,'Torsion')) && any(strcmp(DoF,'InBend')) && ~any(strcmp(DoF,'OutBend')) && ~any(strcmp(DoF,'Axial')) || any(strcmp(DoF,'Torsion')) && ~any(strcmp(DoF,'InBend')) && ~any(strcmp(DoF,'OutBend')) && any(strcmp(DoF,'Axial'))
            C_di = DCM(1,sum(q(n+1:i_element+n)))*C_di;
        elseif any(strcmp(DoF,'Torsion')) && ~any(strcmp(DoF,'InBend')) && ~any(strcmp(DoF,'OutBend')) && ~any(strcmp(DoF,'Axial'))
            C_di = DCM(1,sum(q(1:i_element)))*C_di;
        end

        C_d0 = C_di*C_i0; % Transformation matrix from the global axis to the deformed axis

        if any(strcmp(DoF,'Axial')) && any(strcmp(DoF,'InBend')) && any(strcmp(DoF,'OutBend'))
            r1_d_element = C_d0*r0_0_element + [beam(i_beam).L_element + q(i_element+2*n);0;0];
        elseif any(strcmp(DoF,'Axial')) && ~any(strcmp(DoF,'InBend')) && any(strcmp(DoF,'OutBend')) || any(strcmp(DoF,'InBend')) && ~any(strcmp(DoF,'OutBend'))
            r1_d_element = C_d0*r0_0_element + [beam(i_beam).L_element + q(i_element+n);0;0];
        elseif any(strcmp(DoF,'Axial')) && ~any(strcmp(DoF,'InBend')) && ~any(strcmp(DoF,'OutBend'))
            r1_d_element = C_d0*r0_0_element + [beam(i_beam).L_element + q(i_element);0;0];
        elseif ~any(strcmp(DoF,'Axial'))
            r1_d_element = C_d0*r0_0_element + [beam(i_beam).L_element;0;0];
        end
        
        beam(i_beam).r0(i_element,:) = r0_0_element.';
        beam(i_beam).r1(i_element,:) = (C_d0.'*r1_d_element).';
        beam(i_beam).rCM(i_element,:) = (beam(i_beam).r0(i_element,:) + beam(i_beam).r1(i_element,:))/2;
        r0_0_element = C_d0.'*r1_d_element;

        %         if i_element ~= 1
        %             r1x(i_element) = r0x(i_element) + (beam(i_beam).L_element + q(2*n+i_element))*cos(deg2rad(beam(i_beam).Gamma_deg) + sum(q(1:i_element)))*cos(deg2rad(beam(i_beam).Lambda_deg) + sum(q(n+1:n+i_element)));% + q(i_element-1)); %+ q(i_element))*cos(deg2rad(beam(i_beam).Gamma_deg) + sum(q(2*n+1:2*n+i_element)))*cos(deg2rad(beam(i_beam).Lambda_deg) + sum(q(3*n+1:3*n+i_element)));
        %             r1y(i_element) = r0y(i_element) + (beam(i_beam).L_element + q(2*n+i_element))*cos(deg2rad(beam(i_beam).Gamma_deg) + sum(q(1:i_element)))*sin(deg2rad(beam(i_beam).Lambda_deg) + sum(q(n+1:n+i_element))); %+ (beam(i_beam).L_element + q(i_element))*cos(deg2rad(beam(i_beam).Gamma_deg) + sum(q(2*n+1:2*n+i_element)))*sin(deg2rad(beam(i_beam).Lambda_deg) + sum(q(3*n+1:3*n+i_element)));
        %             r1z(i_element) = r0z(i_element) + (beam(i_beam).L_element + q(2*n+i_element))*sin(deg2rad(beam(i_beam).Gamma_deg) + sum(q(1:i_element)));% + q(i_element-1));% + q(i_element))*sin(deg2rad(beam(i_beam).Gamma_deg) + sum(q(2*n+1:2*n+i_element)));
        %         else
        %             r1x(i_element) = r0x(i_element) + (beam(i_beam).L_element + q(2*n+i_element))*cos(deg2rad(beam(i_beam).Gamma_deg) + q(i_element))*cos(deg2rad(beam(i_beam).Lambda_deg) + q(n+i_element)); %+ q(i_element))*cos(deg2rad(beam(i_beam).Gamma_deg) + sum(q(2*n+1:2*n+i_element)))*cos(deg2rad(beam(i_beam).Lambda_deg) + sum(q(3*n+1:3*n+i_element)));
        %             r1y(i_element) = r0y(i_element) + (beam(i_beam).L_element + q(2*n+i_element))*cos(deg2rad(beam(i_beam).Gamma_deg) + q(i_element))*sin(deg2rad(beam(i_beam).Lambda_deg) + q(n+i_element)); %+ (beam(i_beam).L_element + q(i_element))*cos(deg2rad(beam(i_beam).Gamma_deg) + sum(q(2*n+1:2*n+i_element)))*sin(deg2rad(beam(i_beam).Lambda_deg) + sum(q(3*n+1:3*n+i_element)));
        %             r1z(i_element) = r0z(i_element) + (beam(i_beam).L_element + q(2*n+i_element))*sin(deg2rad(beam(i_beam).Gamma_deg) + q(i_element));% + q(i_element))*sin(deg2rad(beam(i_beam).Gamma_deg) + sum(q(2*n+1:2*n+i_element)));
        %         end
        %         rCMx(i_element) = (r0x(i_element) + r1x(i_element))/2;
        %         rCMy(i_element) = (r0y(i_element) + r1y(i_element))/2;
        %         rCMz(i_element) = (r0z(i_element) + r1z(i_element))/2;
        %
        %         if i_element < beam(i_beam).n_element
        %             r0x(i_element+1) = r1x(i_element);
        %             r0y(i_element+1) = r1y(i_element);
        %             r0z(i_element+1) = r1z(i_element);
        %         else
        %         end
    end
    %     beam(i_beam).r0 = [r0x.' r0y.' r0z.'];
    %     beam(i_beam).r1 = [r1x.' r1y.' r1z.'];
    %     beam(i_beam).rCM = [rCMx.' rCMy.' rCMz.'];
    rCM = reshape(beam(i_beam).rCM,[],1);
end
end