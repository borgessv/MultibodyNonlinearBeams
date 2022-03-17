clear
close all
clc

n_element = 20;
n_frame = 100;
test1 = linspace(0,pi/6,20).';
test2 = linspace(pi/6,pi/4,10).';
test3 = linspace(pi/4,pi/3,10).';
q1 = 0*[0*linspace(0,0.1,20).';-0*test1;0.5*test1;0*test1];
q2 = 0*[0*linspace(0,0.1,20).';-0*test1;0.5*test1;0*test1];
q3 = 0*[0*linspace(0,0.1,10).';-0*test1;0*test1;0*test1];
q4 = 0*[0*linspace(0,0.1,10).';-0*test1;0*test1;0*test1];
q5 = 0*[0*rand(20,1);-2*test1;-test1;0.5*test1];
q6 = 0*[0*rand(10,1);-2*test2;-test2;0.5*test2];
q7 = 0*[0*rand(10,1);-2*test3;-test3;0.5*test3];
q8 = 0*[0*rand(20,1);2*test1;-test1;-0.5*test1];
q9 = 0*[0*rand(10,1);2*test2;-test2;-0.5*test2];
q10 = 0*[0*rand(10,1);2*test3;-test3;-0.5*test3];

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
t = sin(linspace(0,2*pi,n_frame));
v = VideoWriter('video.avi');
v.Quality = 100;
v.FrameRate = 10;
open(v);
f = figure(1);
f.WindowState = 'maximized';
for i = 1:n_frame
    plot_structure_test(Q,'beams.mat','DisplayElements','DisplayUndeformed')
    hold off
    xlim([-25 15])
    ylim([-35 35])
    zlim([-20 20])
    view(45,25)
    %drawnow
    
    frame(i) = getframe(gcf);
    writeVideo(v,frame(i));
   
    q1 = [t(i)*0*linspace(0,0.2,10).';t(i)*0*test2;t(i)*0.25*test2;abs(t(i)*0*test2)];
    q2 = [t(i)*0*linspace(0,0.2,20).';t(i)*0*test1;t(i)*0.25*test1;abs(t(i)*0*test1)];
    q3 = [t(i)*0*linspace(0,0.2,10).';t(i)*0*test2;t(i)*0.25*test3;abs(t(i)*0*test2)];
    q4 = [t(i)*0*linspace(0,0.2,10).';t(i)*0*test2;t(i)*0.25*test2;abs(t(i)*0*test2)];
    q5 = [t(i)*0.5*linspace(0,0.2,20).';t(i)*-2*test1;t(i)*-test1;abs(t(i)*0.5*test1)];
    q6 = [t(i)*0.5*linspace(0,0.2,10).';t(i)*-2*test2;t(i)*-test2;abs(t(i)*0.5*test2)];
    q7 = [t(i)*0.5*linspace(0,0.2,10).';t(i)*-2*test3;t(i)*-test3;abs(t(i)*0.5*test3)];
    q8 = [t(i)*0.5*linspace(0,0.2,20).';t(i)*2*test1;t(i)*-test1;-abs(t(i)*0.5*test1)];
    q9 = [t(i)*0.5*linspace(0,0.2,10).';t(i)*2*test2;t(i)*-test2;-abs(t(i)*0.5*test2)];
    q10 = [t(i)*0.5*linspace(0,0.2,10).';t(i)*2*test3;t(i)*-test3;-abs(t(i)*0.5*test3)];
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
end
close (v)
%%
% fig = figure(2);
% set(gcf, 'Position',  [250, 50, 800, 600])
% movie(fig,frame,10,20)
