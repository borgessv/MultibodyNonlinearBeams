clear
close all
clc

n_element = 20;
n_frame = 200;
test = linspace(pi/36,pi/6,n_element).';
q = 0*[linspace(0,0.2,n_element).';-3*test;-test;0.5*test];
t = sin(0:2*pi/n_frame:2*pi);
v = VideoWriter('video.avi');
v.Quality = 100;
v.FrameRate = 10;
open(v);
f = figure(1);
f.WindowState = 'maximized';
for i = 1:n_frame
    plot_structure(q,'NACA_2412.txt',30,[zeros(3,n_element*2)],linspace(3,2,n_element),-0.5,'DisplayElements','DisplayUndeformed');
    xlim([-0.5 32])
    ylim([-1 6])
    zlim([-10 10])
    view(45,30)
    %drawnow
    
    frame(i) = getframe(gcf);
    writeVideo(v,frame(i));
   
    q = [t(i)*linspace(0,0.2,n_element).';t(i)*-3*test;t(i)*-test;abs(t(i)*0.5*test)];
end
close (v)
%%
% fig = figure(2);
% set(gcf, 'Position',  [250, 50, 800, 600])
% movie(fig,frame,10,20)
