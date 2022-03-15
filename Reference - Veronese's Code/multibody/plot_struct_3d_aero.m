function [Xvec,Zvec_T] = plot_struct_3d_aero(L0,c,a,L_vec,X)
N = length(X)/2;
%Bending
%Torsion
Xvec = zeros(N,1);
Zvec = zeros(N,1);
Zvec_T = zeros(N,1);
Yvec_T = zeros(N,1);
Zvec_T2 = zeros(N,1);
Yvec_T2 = zeros(N,1);
Zvec_T(1,1) = 0;
Yvec_T(1,1) = 0;
Zvec_T2(1,1) = 0;
Yvec_T2(1,1) = 0;
nx_vec = zeros(N,1);
ny_vec = zeros(N,1);
nz_vec = zeros(N,1);
L_vec_x = zeros(N,1);
L_vec_y = zeros(N,1);
L_vec_z = zeros(N,1);
qq = 0;
qq_T = 0;
scale = double(0.8*L0*N/max(L_vec));
for i = 1:N
    qq = qq + X(i,1);
    qq_T = qq_T + X(i+N,1);
    nx_vec(i) = 0;
    ny_vec(i) = 0;
    nz_vec(i) = L_vec(i);
    Xvec(i+1,1) = Xvec(i,1) + L0*cos(qq);
    Zvec(i+1,1) = Zvec(i,1) + L0*sin(qq);
    Yvec_T(i+1,1) = (0.5*c -a)*cos(qq_T);
    Zvec_T(i+1,1) = (0.5*c -a)*sin(qq_T);
    Yvec_T2(i+1,1) = (0.5*c +a)*cos(qq_T);
    Zvec_T2(i+1,1) = (0.5*c +a)*sin(qq_T);
    L_vec_x(i) = Xvec(i)+0.5*L0*cos(qq);
    L_vec_y(i) = ((0.5*c-a)-0.25*c)*cos(qq_T);
    L_vec_z(i) = Zvec(i) +(-(0.5*c-a)+0.25*c)*sin(qq_T);
end
hold on
for i = 1:N
Yvecplot_T = linspace(-Yvec_T2(i+1),Yvec_T(i+1),2);
Xvecplot_T = linspace(Xvec(i),Xvec(i+1),2);
Zvecplot_T = [-Zvec_T2(i+1)+ Zvec(i),Zvec_T(i+1)+ Zvec(i)...
    ;-Zvec_T2(i+1)+ Zvec(i+1),Zvec_T(i+1)+ Zvec(i+1)];
[Xx_T,Yy_T] = meshgrid(Xvecplot_T,Yvecplot_T);
%[Xx_T,Zz_T] = meshgrid(Xvecplot_T,Zvecplot_T);
surf(Xx_T,Yy_T,transpose(Zvecplot_T),'FaceColor','#FFC367');
quiver3(L_vec_x(i),L_vec_y(i),L_vec_z(i),nx_vec(i),ny_vec(i),nz_vec(i),scale,'b')
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
zlabel('z position (m)')
axis('equal')
grid on
xlim([0*N*L0 N*L0])
ylim([-0.25*N*L0 0.25*N*L0])
zlim([-0.3*N*L0 0.85*N*L0])
view([37.5 30]);
end
