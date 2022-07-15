function animation(T,X,H,DoF,element,filename,varargin)
progressbar('loading animation...');

global beam
addpath CrossSectionData

if isnumeric(varargin{end})
    X_r = varargin{end-1};
    H_r = varargin{end};
end

n_frame = length(T);
f = figure('visible','off');
%f.WindowState = 'maximized';
f.Position = [150, 42, 950, 645];
f.Color = 'w';

if any(strcmp(varargin,'gif'))
    gif_filename = [filename '.gif'];
    % % Prellocates variables used to create .gif if MATLAB version is prior to 2022a release
    % frame = struct('cdata',0,'colormap',0);
    % im = cell(n_frame,1);
end
if any(strcmp(varargin,'avi'))
    avi_filename = [filename '.avi'];
    v = VideoWriter(avi_filename);
    v.Quality = 100;
    if T(end) - T(1) < 1
        id = 0.1*length(T);
    else
        [~,id] = min(abs(T - (T(1) + 1)));
    end
    v.FrameRate = id - 1;
    open(v);
    ax = subplot(3,2,[1,3]);
    ax3 = subplot(3,2,2);
    ax2 = subplot(3,2,4);
    ax4 = subplot(3,2,6);
end
if isempty(varargin)
    error('File extension not supported or not specified! Valid options are: "gif", "avi" or both.')
end

r_tip = zeros(n_frame,3);
r_tip_ROM = zeros(n_frame,3);
for i = 1:n_frame
    t1 = tic;
    element_position(X(i,:),DoF);
    for i_beam = 1:length(beam)
        r = [beam(i_beam).r0(1,:);beam(i_beam).r1];
        rCM = beam(i_beam).rCM;
        r_tip(i,:) = r(end,:);
        if any(strcmp(element,'1D'))
            if i == 1
                fig(1) = plot3(ax,r(:,1),r(:,2),r(:,3),'-|k','markersize',2,'markerfacecolor','k');
                hold on
                fig(2) = plot3(ax,rCM(:,1),rCM(:,2),rCM(:,3),'sr','markersize',2,'markerfacecolor','r');
                xlabel('x [m]','interpreter','latex')
                ylabel('y [m]','interpreter','latex')
                zlabel('z [m]','interpreter','latex')
                grid on
                axis equal
                xlim([-6 16])
                ylim([-5 5])
                zlim([-14 14])
            else
                set(0, 'CurrentFigure', f)
                set(fig(1),'XData',r(:,1),'YData',r(:,2),'ZData',r(:,3));
                set(fig(2),'XData',rCM(:,1),'YData',rCM(:,2),'ZData',rCM(:,3));
            end

            if length(varargin) > 2
                element_position(X_r(i,:),DoF);
                r_ROM = [beam(i_beam).r0(1,:);beam(i_beam).r1];
                rCM_ROM = beam(i_beam).rCM;
                if i == 1
                    fig(3) = plot3(ax,r_ROM(:,1),r_ROM(:,2),r_ROM(:,3),'-|b','markersize',2,'markerfacecolor','b');
                    fig(4) = plot3(ax,rCM_ROM(:,1),rCM_ROM(:,2),rCM_ROM(:,3),'sg','markersize',2,'markerfacecolor','g');
                else
                    set(0, 'CurrentFigure', f)
                    set(fig(3),'XData',r_ROM(:,1),'YData',r_ROM(:,2),'ZData',r_ROM(:,3));
                    set(fig(4),'XData',rCM_ROM(:,1),'YData',rCM_ROM(:,2),'ZData',rCM_ROM(:,3));
                end
                legend('FOM','','ROM','','Position',[0.68 0.88 0.05 0.05])
                str = {['Time elapsed: ',num2str(T(i),'%.3f'),' s'],['Tip Position: ','$\hspace{0.7cm}$','FOM','$\hspace{1.8cm}$','ROM'],['$\hspace{1cm}$','x-axis: ',num2str(r(end,1),'%.4e'),' m','$\hspace{0.5cm}$',num2str(r_ROM(end,1),'%.4e'),' m'],['$\hspace{1cm}$','y-axis: ',num2str(r(end,2),'%.4e'),' m','$\hspace{0.5cm}$',num2str(r_ROM(end,2),'%.4e'),' m'],['$\hspace{1cm}$','z-axis: ',num2str(r(end,3),'%.4e'),' m','$\hspace{0.5cm}$',num2str(r_ROM(end,3),'%.4e'),' m']};
                annotation('textbox',[0.2, 0.85, 0.1, 0.1],'String',str,'edgecolor','w','FitBoxToText','on','FontSize',10,'interpreter','latex');
            else
                str = {['Time elapsed: ',num2str(T(i),'%.3f'),' s'],['Tip Position: x-axis: ',num2str(r(end,1),'%.3f'),' m'],['$\hspace{2.15cm}$','y-axis: ',num2str(r(end,2),'%.3f'),' m'],['$\hspace{2.15cm}$','z-axis: ',num2str(r(end,3),'%.3f'),' m']};
                annotation('textbox',[0.2, 0.85, 0.1, 0.1],'String',str,'edgecolor','w','FitBoxToText','on','FontSize',10,'interpreter','latex');
            end


        elseif any(strcmp(element,'3D'))
            Xsection = beam(i_beam).Xsection;
            if endsWith(Xsection,'.txt') == 1
                Xsection_data = readtable(Xsection);
                Xsection_data = Xsection_data{:,:};
                Xsection_data(:,3) = zeros(length(Xsection_data),1);
            elseif any(strcmp(Xsection,'Circular'))
                angle_circular = 0:pi/90:2*pi;
                x_circular = cos(angle_circular)/2 + 0.5;
                y_circular = sin(angle_circular)/2;
                Xsection_data = [x_circular.' y_circular.' zeros(length(x_circular),1)];
            elseif any(strcmp(Xsection,'Rectangular'))
                Xsection_data = [0 -beam(i_beam).t_beam/2; 0 beam(i_beam).t_beam/2; 1 beam(i_beam).t_beam/2; 1 -beam(i_beam).t_beam/2; 0 -beam(i_beam).t_beam/2];
                Xsection_data(:,3) = zeros(length(Xsection_data),1);
            end
            for i_element = 1:beam(i_beam).n_element
                Xsection_0 = beam(i_beam).element(i_element).C_d0'*[Xsection_data(:,3) beam(i_beam).element(i_element).c*(Xsection_data(:,1)-beam(i_beam).yCM) beam(i_beam).element(i_element).c*(Xsection_data(:,2)-beam(i_beam).zCM)].';
                r0_Xsection_0_def = (r(i_element,:) + Xsection_0.').';
                r1_Xsection_0_def = (r(i_element+1,:) + Xsection_0.').';
                X3D_def(:,2*i_element-1:2*i_element) = [r0_Xsection_0_def(1,:).' r1_Xsection_0_def(1,:).'];
                Y3D_def(:,2*i_element-1:2*i_element) = [r0_Xsection_0_def(2,:).' r1_Xsection_0_def(2,:).'];
                Z3D_def(:,2*i_element-1:2*i_element) = [r0_Xsection_0_def(3,:).' r1_Xsection_0_def(3,:).'];
            end
            if i == 1
                % Instant Position:
                fig1 = surf(ax,X3D_def,Y3D_def,Z3D_def,'facecolor',[1 0.5 0],'linestyle','none');
                hold(ax,'on')
                fig2 = plot3(ax,X3D_def(:),Y3D_def(:),Z3D_def(:),'color',[0.5 0.5 0.5]);
                title(ax,'Structure Simulation',['Time Elapsed: ',num2str(T(i),'%.3f'),' s'],'interpreter','latex','fontsize',10)
                xlabel(ax,'x [m]','interpreter','latex','FontSize',9)
                ylabel(ax,'y [m]','interpreter','latex','FontSize',9)
                zlabel(ax,'z [m]','interpreter','latex','FontSize',9)
                grid(ax,'on')
                axis(ax,'equal')
                xlim(ax,[-6 16])
                ylim(ax,[-12 12])
                zlim(ax,[-14 14])
                %view(ax,0,25)
                
                % Tip Position:
                fig11 = plot(ax3,T(1),r_tip(1,1),'-b','linewidth',1.5);
                hold(ax3,'on')
                fig12 = plot(ax3,T(1),r_tip(1,2),'-r','linewidth',1.5);
                fig13 = plot(ax3,T(1),r_tip(1,3),'-k','linewidth',1.5);
                %title(ax3,'Tip Position','interpreter','latex','FontSize',9)
                xlabel(ax3,'t [s]','interpreter','latex','FontSize',9)
                ylabel(ax3,'Tip Position [m]','interpreter','latex','FontSize',9)
                grid(ax3,'on')
                xlim(ax3,[0 T(end)])
                ylim(ax3,[-16 16])

                % Tip Twist Angle:
                if any(strcmp(DoF,'Torsion'))
                    fig14 = plot(ax2,T(1),rad2deg(sum(X(1,end-beam(i_beam).n_element + 1:end))),'-b','linewidth',1.5);
                    hold(ax2,'on')
                else
                    fig14 = plot(ax2,T(1),0,'-b','linewidth',1.5);
                    hold(ax2,'on')
                end
                xlabel(ax2,'t [s]','interpreter','latex','FontSize',9)
                ylabel(ax2,'Tip Twist Angle [deg]','interpreter','latex','FontSize',9)
                grid(ax2,'on')
                xlim(ax2,[0 T(end)])
                %ylim(ax2,[-16 16])
                
                % Hamiltonian:
                fig15 = plot(ax4,T(1),H(1),'-r','linewidth',1.5);
                hold(ax4,'on')
                xlabel(ax4,'t [s]','interpreter','latex','FontSize',9)
                ylabel(ax4,'$\mathcal{H}$ [J]','interpreter','latex','FontSize',9)
                grid(ax4,'on')
                xlim(ax4,[0 T(end)])
                %ylim(ax4,[-16 16])

            else
                set(0, 'CurrentFigure', f)
                set(fig1,'XData',X3D_def,'YData',Y3D_def,'ZData',Z3D_def);
                title(ax,'Structure Simulation',['Time Elapsed: ',num2str(T(i),'%.3f'),' s'],'interpreter','latex')
                set(fig2,'XData',X3D_def(:),'YData',Y3D_def(:),'ZData',Z3D_def(:));
                set(fig11,'XData',T(1:i),'YData',r_tip(1:i,1));
                set(fig12,'XData',T(1:i),'YData',r_tip(1:i,2));
                set(fig13,'XData',T(1:i),'YData',r_tip(1:i,3));
                if any(strcmp(DoF,'Torsion'))
                    set(fig14,'XData',T(1:i),'YData',rad2deg(sum(X(1:i,end - beam(i_beam).n_element + 1:end),2)));
                else
                    set(fig14,'XData',T(1:i),'YData',zeros(1,i));
                end
                set(fig15,'XData',T(1:i),'YData',H(1:i));
                %view(ax,360*(i-1)/n_frame,25)
            end
            if isnumeric(varargin{end})
                element_position(X_r(i,:),DoF);
                r_ROM = [beam(i_beam).r0(1,:);beam(i_beam).r1];
                r_tip_ROM(i,:) = r_ROM(end,:);
                for i_element = 1:beam(i_beam).n_element
                    Xsection_0 = beam(i_beam).element(i_element).C_d0'*[Xsection_data(:,3) beam(i_beam).element(i_element).c*(Xsection_data(:,1)-beam(i_beam).yCM) beam(i_beam).element(i_element).c*(Xsection_data(:,2)-beam(i_beam).zCM)].';
                    r0_Xsection_0_def = (r_ROM(i_element,:) + Xsection_0.').';
                    r1_Xsection_0_def = (r_ROM(i_element+1,:) + Xsection_0.').';
                    X3D_def(:,2*i_element-1:2*i_element) = [r0_Xsection_0_def(1,:).' r1_Xsection_0_def(1,:).'];
                    Y3D_def(:,2*i_element-1:2*i_element) = [r0_Xsection_0_def(2,:).' r1_Xsection_0_def(2,:).'];
                    Z3D_def(:,2*i_element-1:2*i_element) = [r0_Xsection_0_def(3,:).' r1_Xsection_0_def(3,:).'];
                end
                if i == 1
                    % Structure's instant position:
                    fig3 = surf(ax,X3D_def,Y3D_def,Z3D_def,'facecolor',[0.3010 0.7450 0.9330],'linestyle','none');
                    fig4 = plot3(ax,X3D_def(:),Y3D_def(:),Z3D_def(:),'color',[0.5 0.5 0.5]);
                    legend(ax,[fig1 fig3],'FOM','ROM','location','northeast','interpreter','latex','FontSize',9);%,[0.68 0.88 0.05 0.05])
                    
                    % % Equilibrium Position:
                    % surf(ax2,X3D_def,Y3D_def,Z3D_def,'facecolor',[0.3010 0.7450 0.9330],'linestyle','none');
                    % plot3(ax2,X3D_def(:),Y3D_def(:),Z3D_def(:),'color',[0.5 0.5 0.5]);
                    
                    % Tip Position:
                    fig21 = plot(ax3,T(1),r_ROM(end,1),'-c','linewidth',1.5);
                    fig22 = plot(ax3,T(1),r_ROM(end,2),'-m','linewidth',1.5);
                    fig23 = plot(ax3,T(1),r_ROM(end,3),'-','color',[0.5 0.5 0.5],'linewidth',1.5);
                    legend([fig11 fig12 fig13 fig21 fig22 fig23],'x - FOM','y - FOM','z - FOM','x - ROM','y - ROM','z - ROM','location','northeast','interpreter','latex','FontSize',8);

                    % Tip Twist Angle:
                    if any(strcmp(DoF,'Torsion'))
                        fig24 = plot(ax2,T(1),rad2deg(sum(X_r(1,end-beam(i_beam).n_element + 1:end))),'-c','linewidth',1.5);
                    else
                        fig24 = plot(ax2,T(1),0,'-c','linewidth',1.5);
                    end
                    legend([fig14 fig24],'FOM','ROM','location','northeast','interpreter','latex','FontSize',8);

                    % Hamiltonian:
                    fig25 = plot(ax4,T(1),H_r(1),'-m','linewidth',1.5);
                    legend([fig15 fig25],'FOM','ROM','location','northeast','interpreter','latex','FontSize',8);

                else
                    set(0, 'CurrentFigure', f)
                    set(fig3,'XData',X3D_def,'YData',Y3D_def,'ZData',Z3D_def);
                    set(fig4,'XData',X3D_def(:),'YData',Y3D_def(:),'ZData',Z3D_def(:));
                    set(fig21,'XData',T(1:i),'YData',r_tip_ROM(1:i,1));
                    set(fig22,'XData',T(1:i),'YData',r_tip_ROM(1:i,2));
                    set(fig23,'XData',T(1:i),'YData',r_tip_ROM(1:i,3));
                    if any(strcmp(DoF,'Torsion'))
                        set(fig24,'XData',T(1:i),'YData',rad2deg(sum(X_r(1:i,end - beam(i_beam).n_element + 1:end),2)));
                    else
                        set(fig24,'XData',T(1:i),'YData',zeros(1,i));
                    end
                    set(fig25,'XData',T(1:i),'YData',H_r(1:i));
                end
                %str = {['Time elapsed: ',num2str(T(i),'%.3f'),' s']};%,['Tip Position: ','$\hspace{0.7cm}$','FOM','$\hspace{1.8cm}$','ROM'],['$\hspace{1cm}$','x-axis: ',num2str(r(end,1),'%.4e'),' m','$\hspace{0.5cm}$',num2str(r_ROM(end,1),'%.4e'),' m'],['$\hspace{1cm}$','y-axis: ',num2str(r(end,2),'%.4e'),' m','$\hspace{0.5cm}$',num2str(r_ROM(end,2),'%.4e'),' m'],['$\hspace{1cm}$','z-axis: ',num2str(r(end,3),'%.4e'),' m','$\hspace{0.5cm}$',num2str(r_ROM(end,3),'%.4e'),' m']};
                %annotation('textbox',[0.1, 0.85, 0.1, 0.1],'String',str,'LineStyle','none','FitBoxToText','on','FontSize',9,'interpreter','latex');
            else
                %str = {['Time elapsed: ',num2str(T(i),'%.3f'),' s']};%,['Tip Position: x-axis: ',num2str(r(end,1),'%.3f'),' m'],['$\hspace{2.15cm}$','y-axis: ',num2str(r(end,2),'%.3f'),' m'],['$\hspace{2.15cm}$','z-axis: ',num2str(r(end,3),'%.3f'),' m']};
                %annotation('textbox',[0.1, 0.85, 0.1, 0.1],'String',str,'edgecolor','w','FitBoxToText','on','FontSize',9,'interpreter','latex');
                legend([fig11 fig12 fig13],'x','y','z','location','northeast','interpreter','latex','fontsize',8);
            end
        end
    end
    

    if any(strcmp(varargin,'gif'))
        if i == 1
            exportgraphics(f,gif_filename,'Append',false)
        else
            exportgraphics(f,gif_filename,'Append',true)
        end
        % % Creates .gif for MATLAB versions prior to 2022a:
        % frame(i) = getframe(gcf); % Uncomment for .gif if MATLAB version is prior to 2022a release
        % im{i} = frame2im(frame(i));
        % [A,map] = rgb2ind(im{i},256);
        % if i == 1
        %    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0);
        % elseif any(i == 2:n_frame)
        %    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0);
        % end
    end
    if any(strcmp(varargin,'avi'))
        writeVideo(v,getframe(f));
    end

    % Updates progress bar:
    t2 = toc(t1);
    t_rem = t2*(length(T) - i);
    if i == n_frame
        progressbar(i/length(T)*100);
    else
        progressbar(i/length(T)*100,t_rem);
    end
    %delete(findall(f,'type','annotation')) % clears annotation to update in the next iteration
end
if any(strcmp(varargin,'gif'))
    winopen(gif_filename);
end
if any(strcmp(varargin,'avi'))
    close(v)
end
end