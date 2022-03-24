clear
close all
clc

n_frame = 100;
test1 = linspace(0,pi/6,20).';
test2 = linspace(pi/6,pi/4,10).';
test3 = linspace(pi/4,pi/3,10).';
test4 = linspace(0,pi/6,10).';
t = sin(linspace(0,2*pi,n_frame));

% % Video settings:
% v = VideoWriter('video.avi');
% v.Quality = 100;
% v.FrameRate = 50;
% open(v);

% Figure settings:
f = figure('visible','off');
% f.WindowState = 'maximized';
set(gcf, 'Position',  [250, 42, 750, 645])
set(gcf,'color','w');

% Gif settings:
filename = 'simulation.gif';

progressbar('loading animation...');
frame = struct('cdata',0,'colormap',0);
im = cell(n_frame,1);
for i = 1:n_frame
%     t1 = tic;
    q1 = [t(i)*0*linspace(0,0.2,10).';t(i)*0*test2;t(i)*0.25*test2;abs(t(i)*0*test2)];
    q2 = [t(i)*0*linspace(0,0.2,20).';t(i)*0*test1;t(i)*0.25*test1;abs(t(i)*0*test1)];
    q3 = [t(i)*0*linspace(0,0.2,10).';t(i)*0*test2;t(i)*0.25*test3;abs(t(i)*0*test2)];
    q4 = [t(i)*0*linspace(0,0.2,10).';t(i)*0*test2;t(i)*0.25*test2;abs(t(i)*0*test2)];
    q5 = [t(i)*linspace(0,0.3,20).';-2*t(i)*test1;t(i)*-test1;abs(t(i)*0.5*test1)];
    q6 = [t(i)*linspace(0,0.3,10).';-2*t(i)*test2;t(i)*-test2;abs(t(i)*0.5*test2)];
    q7 = [t(i)*linspace(0,0.3,10).';-2*t(i)*test3;t(i)*-test3;abs(t(i)*0.5*test3)];
    q8 = [t(i)*linspace(0,0.3,20).';2*t(i)*test1;t(i)*-test1;-abs(t(i)*0.5*test1)];
    q9 = [t(i)*linspace(0,0.3,10).';2*t(i)*test2;t(i)*-test2;-abs(t(i)*0.5*test2)];
    q10 = [t(i)*linspace(0,0.3,10).';2*t(i)*test3;t(i)*-test3;-abs(t(i)*0.5*test3)];
    q11 = [t(i)*linspace(0,0.2,10).';-t(i)*test4;-t(i)*test4;abs(t(i)*0.5*test4)];
    q12 = [t(i)*linspace(0,0.2,10).';t(i)*test4;-t(i)*test4;-abs(t(i)*0.5*test4)];
    q13 = [t(i)*linspace(0,0.1,10).';t(i)*0*test4;abs(0*t(i)*test4);t(i)*test4];
    Q(1).q = q1;
    Q(2).q = q2;
    Q(3).q = q3;
    Q(4).q = q4;
    Q(5).q = q5;
    Q(6).q = q6;
    Q(7).q = q7;
    Q(8).q = q8;
    Q(9).q = q9;
    Q(10).q = q10;
    Q(11).q = q11;
    Q(12).q = q12;
    Q(13).q = q13;
    
    plot_structure_test(Q,'beams.mat','DisplayElements','DisplayUndeformed')
    hold off
    xlim([-25 15])
    ylim([-35 35])
    zlim([-25 25])
    view(60,20)
    frame(i) = getframe(gcf);

%     % Creating .avi video:
%     writeVideo(v,frame(i));

    % Creating .gif:
    im{i} = frame2im(frame(i));
    [A,map] = rgb2ind(im{i},256);
    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0);
    elseif any(i == 0:4:n_frame)
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0);
    end

    progressbar(i/n_frame*100);
    %waitbar(i/n_frame)
%     t2 = toc(t1);
%     t_rem = n_frame*(t2) - i*(t2);
%     fprintf('Remaining time: %.2f',t_rem)
end
% close(v)
web('simulation.gif')
