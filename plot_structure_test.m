function plot_structure_test(Q,beam_data,Xsection,varargin)
%
% plot_structure(q,Xsection,L,r0_ea_v,c,a,varargin) plots the deformed 
% beam in a multibody analysis considering axial, in-plane bending, 
% out-of-plane bending and torsion displacements.
% 
% Inputs:
%       - q: vector of generalized coordinates in the sequence [q_axial;
%       q_torsion; q_outofplane; q_inplane] such that the size of q is
%       4*number-of-elements;
%       - Xsection: coordinate data file (.txt) of the cross-section
%       considered in the analysis. The inputs 'Rectangular' and 'Circular'
%       can be used instead, but if 'Rectangular' the tickness of the first
%       elementmust be informed in varargin{1};
%       - L: undeformed total length/span of the beam;
%       - r0_ea_v: initial/undeformed position of each element's elastic 
%       axis. Must have the same size as q; 
%       - c: chord/width/diameter of the beam. Must be specified for each 
%       element, therefore its size must be equal the number of elements;
%       - a: elastic axis position along the chord/width of the beam. Must
%       be in the range [-1,1], where 0 is the center, i.e. c/2;
%       - varargin: optional arguments:
%           - d: tickness of the first element if 'Rectangular' is used 
%           for Xsection variable. 
%           Must be the first argument of varargin;
%           - 'DisplayElements': plot the contour of each element;
%           - 'DisplayUndeformed': plot the undeformed beam;
% 
% Outputs:
%       - 3D plot of the deformed and undeformed (optional) beam

load(beam_data,'beam')

% if a > 1 || a < -1
%     error('The elastic axis position along the chord, i.e. "a", must be declared in the range [-1,1]')
% else
% end
plot_options = varargin(cellfun(@(x_plt) ischar(x_plt),varargin));

switch isnumeric(varargin{1})
    case 0
        if endsWith(Xsection,'.txt') == 1
            Xsection_data = readtable(Xsection);
            Xsection_data = Xsection_data{:,:};
        elseif any(strcmp(Xsection,'Circular'))
            angle_circular = 0:pi/18:2*pi;
            x_circular = cos(angle_circular)/2 + 0.5;
            y_circular = sin(angle_circular)/2;
            Xsection_data = [x_circular.' y_circular.'];
        else
            error(['Cross-section .txt data file is missing. For ''Rectangular''' ...
                ' cross-section the tickness of the first element must be specified in varargin{1}.'])
        end
    case 1
        d = varargin{1};
        if any(strcmp(Xsection,'Rectangular'))
            Xsection_data = [0 -d/2; 0 d/2; 1 d/2; 1 -d/2; 0 -d/2];
        else
            error(['varargin{1} is numeric. You must delete it in order to use ' ...
                'a .txt data file or specify ''Rectangular'' for the cross-section input variable'])
        end
end

for i_beam = 1:length(beam)

    for i_element = 1:beam(i_beam).n_element
        r0_i = beam(i_beam).C_i0*beam(i_beam).element(i_element).r0;
        r1_i = beam(i_beam).C_i0*beam(i_beam).element(i_element).r1;

        if any(strcmp(plot_options,'DisplayUndeformed')) || any(strcmp(plot_options,'displayundeformed')) || any(strcmp(plot_options,'DISPLAYUNDEFORMED'))
        Xsection_coord_undeformed = [zeros(1,length(Xsection_data));(beam(i_beam).element(i_element).c.*Xsection_data(:,1)-beam(i_beam).element(i_element).c*beam(i_beam).yCM).';(beam(i_beam).element(i_element).c*Xsection_data(:,2)).'];
        r0_Xsection_0 = beam(i_beam).C_i0.'*(r0_i + Xsection_coord_undeformed);
        r1_Xsection_0 = beam(i_beam).C_i0.'*(r1_i + Xsection_coord_undeformed);
        
        X3D_undeformed = [r0_Xsection_0(1,:).' r1_Xsection_0(1,:).'];
        Y3D_undeformed = [r0_Xsection_0(2,:).' r1_Xsection_0(2,:).'];
        Z3D_undeformed = [r0_Xsection_0(3,:).' r1_Xsection_0(3,:).'];
        surf(X3D_undeformed, Y3D_undeformed, Z3D_undeformed,'facecolor','b','linestyle','none');
        hold on
        if any(strcmp(plot_options,'DisplayElements')) || any(strcmp(plot_options,'displayelements')) || any(strcmp(plot_options,'DISPLAYELEMENTS'))
            plot3(X3D_undeformed,Y3D_undeformed,Z3D_undeformed,'k')
        else
            if i_element == 1
                plot3(X3D_undeformed(:,1),Y3D_undeformed(:,1),Z3D_undeformed(:,1),'k')
            elseif i_element == beam(i_beam).n_element
                plot3(X3D_undeformed(:,2),Y3D_undeformed(:,2),Z3D_undeformed(:,2),'k')
            end
        end
        else
        end
    end 
end
axis equal
xlabel('x_0 [m]')
ylabel('y_0 [m]')
zlabel('z_0 [m]')
end