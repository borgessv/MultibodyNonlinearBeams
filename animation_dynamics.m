function animation_dynamics(T,X,DoF,filename,varargin)

global beam
n_DoF = length(DoF)*sum(vertcat(beam.n_element));
n_frame = length(T);
f = figure('visible','off');
set(gcf, 'Position',  [250, 42, 750, 645])
set(gcf,'color','w');

if any(strcmp(varargin,'gif'))
    gif_filename = [filename '.gif'];
    % f = figure('visible','off');
    % set(gcf, 'Position',  [250, 42, 750, 645])
    % set(gcf,'color','w');

    % % Prellocates variables used to create .avi (or .gif if prior to 2022a release)
    % frame = struct('cdata',0,'colormap',0);
    % im = cell(n_frame,1);
end
if any(strcmp(varargin,'avi'))
    % Video settings:
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
    % Figure settings (optional):
    % f = figure('visible','off');
    % f.WindowState = 'maximized';
    % set(gcf, 'Position',  [250, 42, 750, 645])
    % set(gcf,'color','w');
    frame = struct('cdata',0,'colormap',0);
end
if isempty(varargin)
    error('File extension not supported or not specified! Valid options are: "gif", "avi" or both.')
end

progressbar('loading animation...');
for i = 1:n_frame
    %if any(T(i) == 0:0.05:T(end))
    t1 = tic;
    f.Visible = 'off';
    %set(gcf, 'Position',  [250, 42, 750, 645])
    %set(gcf,'color','w');
    element_position(X(i,n_DoF+1:end),DoF);
    r = [beam(1).r0(1,:);beam(1).r1];
    rCM = beam(1).rCM;
    plot3(r(:,1),r(:,2),r(:,3),'-|k','markersize',2,'markerfacecolor','k')
    hold on
    plot3(rCM(:,1),rCM(:,2),rCM(:,3),'sr','markersize',2,'markerfacecolor','r')
    xlabel('x [m]','interpreter','latex')
    ylabel('y [m]','interpreter','latex')
    zlabel('z [m]','interpreter','latex')
    grid on
    axis equal
    xlim([-6 16])
    ylim([-5 5])
    zlim([-14 14])
    str = {['Elapsed time: ',num2str(round(T(i),3)),' s'],['Tip Position: x-axis: ',num2str(round(r(end,1),3)),' m'],['$\hspace{2.15cm}$','y-axis: ',num2str(round(r(end,2),3)),' m'],['$\hspace{2.15cm}$','z-axis: ',num2str(round(r(end,3),3)),' m']};
    annotation('textbox',[0.2, 0.85, 0.1, 0.1],'String',str,'edgecolor','w','FitBoxToText','on','FontSize',10,'interpreter','latex');
    %view(85,2)
    hold off

    if any(strcmp(varargin,'gif'))
        % Creates .gif for MATLAB versions from 2022a:
        if i == 1
            exportgraphics(gca,gif_filename,'Append',false)
        else
            exportgraphics(gca,gif_filename,'Append',true)
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
        frame(i) = getframe(gcf);
        % Appends new frame in .avi file:
        writeVideo(v,frame(i));
    end

    % Updates progress bar:
    t2 = toc(t1);
    t_rem = t2*(length(T) - i);
    if i == n_frame
        progressbar(i/length(T)*100);
    else
        progressbar(i/length(T)*100,t_rem);
    end
    clf % clears figure
end
if any(strcmp(varargin,'gif'))
    web(gif_filename)
end
if any(strcmp(varargin,'avi'))
    close(v) % Uncomment for .avi creation
end
end