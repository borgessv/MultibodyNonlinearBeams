function beam = create_beam(beam_data)

addpath CrossSectionData

% Reading input file and creating arrays for the properties of the beams:
[beam_properties,Xsection] = xlsread(beam_data);
ID_beam = beam_properties(1,:);
n_element_beam = beam_properties(3,:);
L_beam = beam_properties(4,:);
m0_beam = beam_properties(5,:);
m1_beam = beam_properties(6,:);
c0_beam = beam_properties(7,:);
c1_beam = beam_properties(8,:);
xCM_beam = beam_properties(9,:);
yCM_beam = beam_properties(10,:);
zCM_beam = beam_properties(11,:);
Lambda_beam_deg = beam_properties(12,:);
Gamma_beam_deg = beam_properties(13,:);
AoI_beam_deg = beam_properties(14,:);
connectivity_beam = beam_properties(15,:);
connection_point_beam = beam_properties(16,:);
E_beam = beam_properties(17,:);
G_beam = beam_properties(19,:);

n_beam = length(ID_beam);
beam(1:n_beam) = struct('Xsection','','n_element',0,'L',0,'L_element',0,...
    'C_i0',zeros(3),'yCM',0,'zCM',0,...
    'Lambda_deg',0,'Gamma_deg',0,'AoI_deg',0,'m',0,'c',0,...
    'connectivity',0,'connection_point',0,...
    'connection_element',0,'connection_element_point',0,'E',0,'G',0,...
    'element',struct('r0',zeros(3,1),'r1',zeros(3,1),...
    'Xsection_Iyy',0,'Xsection_Izz',0,'Xsection_Iyz',0,...
    'Xsection_J',0,'Xsection_A',0,...
    'x0_element',0,'x1_element',0,'m',0,'c',0));

syms x
for i_beam = 1:n_beam
    beam(i_beam).Xsection = Xsection{3,i_beam+1};
    Xsection_data_raw = readtable(beam(i_beam).Xsection);
    beam(i_beam).n_element = n_element_beam(i_beam);
    beam(i_beam).L = L_beam(i_beam);
    beam(i_beam).L_element = beam(i_beam).L/beam(i_beam).n_element;
    beam(i_beam).Lambda_deg = Lambda_beam_deg(i_beam);
    beam(i_beam).Gamma_deg = Gamma_beam_deg(i_beam);
    beam(i_beam).AoI_deg = AoI_beam_deg(i_beam);
    beam(i_beam).xCM = xCM_beam(i_beam);
    beam(i_beam).yCM = yCM_beam(i_beam);
    beam(i_beam).zCM = zCM_beam(i_beam);
    beam(i_beam).C_i0 = DCM(1,deg2rad(beam(i_beam).AoI_deg))*DCM(2,deg2rad(beam(i_beam).Gamma_deg))*DCM(3,deg2rad(beam(i_beam).Lambda_deg));
    beam(i_beam).m = m0_beam(i_beam) + m1_beam(i_beam)*x;
    beam(i_beam).c = c0_beam(i_beam) + c1_beam(i_beam)*x;
    beam(i_beam).connectivity = connectivity_beam(i_beam);
    beam(i_beam).connection_point = connection_point_beam(i_beam);
    beam(i_beam).E = E_beam(i_beam);
    beam(i_beam).G = G_beam(i_beam);

    for i_element = 1:beam(i_beam).n_element
        if i_element == 1 && isnan(beam(i_beam).connectivity) && isnan(beam(i_beam).connection_point)
            beam(i_beam).element(i_element).r0 = [0;0;0];
            beam(i_beam).element(i_element).x0_element = 0;
            beam(i_beam).element(i_element).x1_element = beam(i_beam).element(i_element).x0_element + beam(i_beam).L_element;
        elseif i_element == 1 && isnan(beam(i_beam).connectivity) && ~isnan(beam(i_beam).connection_point)
            beam(i_beam).element(i_element).r0 = beam(i_beam).C_i0.'*[beam(i_beam).connection_point;0;0];
            beam(i_beam).element(i_element).x0_element = 0;
            beam(i_beam).element(i_element).x1_element = beam(i_beam).element(i_element).x0_element + beam(i_beam).L_element;
        elseif i_element == 1 && ~isnan(beam(i_beam).connectivity)
            beam(i_beam).element(i_element).r0 = beam(beam(i_beam).connectivity).C_i0.'*(beam(beam(i_beam).connectivity).C_i0*beam(beam(i_beam).connectivity).element(1).r0 + [beam(i_beam).connection_point;0;0]);
            beam(i_beam).element(i_element).x0_element = 0;
            beam(i_beam).element(i_element).x1_element = beam(i_beam).element(i_element).x0_element + beam(i_beam).L_element;
        else
            beam(i_beam).element(i_element).r0 = beam(i_beam).element(i_element-1).r1;
            beam(i_beam).element(i_element).x1_element = beam(i_beam).element(i_element).x0_element + beam(i_beam).L_element;

        end
        if i_beam ~= 1
            for i = 1:beam(beam(i_beam).connectivity).n_element
                if beam(i_beam).connection_point >= beam(beam(i_beam).connectivity).element(i).x0_element && beam(i_beam).connection_point <= beam(beam(i_beam).connectivity).element(i).x1_element
                    beam(i_beam).connection_element = i;
                    beam(i_beam).connection_element_point = (beam(i_beam).connection_point - beam(beam(i_beam).connectivity).element(i).x0_element)/(beam(beam(i_beam).connectivity).element(i).x1_element - beam(beam(i_beam).connectivity).element(i).x0_element);
                else
                end
            end
        else
        end
        beam(i_beam).element(i_element).r1 = beam(i_beam).C_i0.'*(beam(i_beam).C_i0*beam(i_beam).element(i_element).r0 + [beam(i_beam).L_element;0;0]);

        beam(i_beam).element(i_element).m = double(vpaintegral(beam(i_beam).m,beam(i_beam).element(i_element).x0_element,beam(i_beam).element(i_element).x1_element));
        beam(i_beam).element(i_element).c = double(mean([subs(beam(i_beam).c,x,beam(i_beam).element(i_element).x1_element),subs(beam(i_beam).c,x,beam(i_beam).element(i_element).x0_element)]));

        % Calculation of cross-section properties for each element:
        Xsection_data = (beam(i_beam).element(i_element).c).*Xsection_data_raw{:,:};
        Xsection_data_interp = [interpft(Xsection_data(:,1),500) interpft(Xsection_data(:,2),500)];
        [A,I,~] = polygeom(Xsection_data_interp(:,1),Xsection_data_interp(:,2));
        beam(i_beam).element(i_element).Xsection_A = A(1);
        beam(i_beam).element(i_element).Xsection_Iyy = I(4);
        beam(i_beam).element(i_element).Xsection_Izz = I(5);
        beam(i_beam).element(i_element).Xsection_Iyz = I(6);

        % Torsion constant approximation for solid cross-section - see Kollbrunner, C. F., et al. Torsion in Structures: an engineering approach. 1969.
        beam(i_beam).element(i_element).Xsection_J = beam(i_beam).element(i_element).Xsection_A^4/(40*(beam(i_beam).element(i_element).Xsection_Iyy + beam(i_beam).element(i_element).Xsection_Izz));
        
        beam(i_beam).element(i_element).EIyy = beam(i_beam).E*beam(i_beam).element(i_element).Xsection_Iyy;
        beam(i_beam).element(i_element).EIzz = beam(i_beam).E*beam(i_beam).element(i_element).Xsection_Izz;
        beam(i_beam).element(i_element).GJ = beam(i_beam).G*beam(i_beam).element(i_element).Xsection_J;
        beam(i_beam).element(i_element).EA = beam(i_beam).E*beam(i_beam).element(i_element).Xsection_A;

        if i_element ~= beam(i_beam).n_element
            beam(i_beam).element(i_element+1).x0_element = beam(i_beam).element(i_element).x1_element;
        else
        end
    end
    save beams_test.mat beam
end
