addpath test_results

%% Test 1: Concentrated tip force:
F = 0:15:150;
z = [0 0.0989948 0.192374 0.276274 0.349142 0.411195 0.463578 0.507736 0.545075 0.576823 0.604002];
x = [1 0.994113 0.977556 0.953051 0.923726 0.892222 0.860386 0.829358 0.799759 0.771879 0.745806];

z_true = [0 0.09895287958115183 0.19240837696335075 0.2757853403141361 0.3490837696335078 0.40955497382198947 0.46269633507853397 0.506675392670157 0.5433246073298429 0.5772251308900523 0.6037958115183245];
x_true = [1 0.9933766233766234 0.9766233766233766 0.9520779220779221 0.9232467532467532 0.8916883116883116 0.8601298701298701 0.8289610389610389 0.8001298701298701 0.7724675324675324 0.7463636363636363];


figure
plot(F,z,'-xr','linewidth',1.5,'markersize',12)
hold on
plot(F,z_true,'--ok','linewidth',1.5)
xlabel('$F_{t1,z}$ [N]','FontSize',16,'Interpreter','latex')
ylabel('$z_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 700, 550])

figure
plot(F,x,'-xr','linewidth',1.5,'markersize',12)
hold on
plot(F,x_true,'--ok','linewidth',1.5)
xlabel('$F_{t1,z}$ [N]','FontSize',16,'Interpreter','latex')
ylabel('$x_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 700, 550])


L2 = sqrt(sum((z-z_true).^2));


%% Test 2: Concentrated Tip Moment
M = 0:15:120;

x = [1 0.985178 0.9415 0.871282 0.77822 0.667158 0.543789 0.414302 0.285003];
z = [0 -0.148394 -0.290214 -0.419274 -0.530138 -0.618429 -0.681074 -0.716465 -0.724515];


z_true = [0
    -0.149484536082474
    -0.292783505154639
    -0.422680412371134
    -0.534020618556701
    -0.622680412371134
    -0.683505154639175
    -0.7185567010309278
    -0.727835051546392
    ];
x_true = [1
    0.985510996119017
    0.943078913324709
    0.872703751617076
    0.779560155239327
    0.666752910737387
    0.542561448900388
    0.413195342820181
    0.282794307891332
    ];


figure
plot(M,x,'-xr','linewidth',1.5,'markersize',12)
hold on
plot(M,x_true,'--ok','linewidth',1.5)
xlabel('$M_{t2,y}$ [Nm]','FontSize',16,'Interpreter','latex')
ylabel('$x_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 700, 550])

figure
plot(M,z,'-xr','linewidth',1.5,'markersize',12)
hold on
plot(M,z_true,'--ok','linewidth',1.5)
xlabel('$M_{t2,y}$ [Nm]','FontSize',16,'Interpreter','latex')
ylabel('$z_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 700, 550])

L2 = sqrt(sum((z-z_true.').^2));


%% Beam Model 1 dynamic response:
global beam
for i = 1:size(q_ROM,1)
    element_positionCM(q_ROM(i,:),DoF);
    r = [beam(1).r0(1,:);beam(1).r1];
    rCM = beam(1).rCM;
    r(:,2) = [rCM(1,2);rCM(:,2)];
    r_tip(i,:) = r(end,:);
    ang_tip(i) = rad2deg(sum(q_ROM(i,57:end)));
end

global beam
Xeq = phi_r*eta_eq;
element_positionCM(Xeq,DoF);
r = [beam(1).r0(1,:);beam(1).r1];
rCM = beam(1).rCM;
r(:,2) = [rCM(1,2);rCM(:,2)];
r_tip(i,:) = r(end,:);
ang_tip(i) = sum(Xeq(57:end));


figure
load test3_cost_n10_qoqa.mat tspan r_tip
r1 = r_tip;
plot(tspan,r1(:,3),'--k','linewidth',1.5)
hold on
load test3_cost_r2_qoqa.mat tspan r_tip
r2 = r_tip;
plot(tspan,r2(:,3),'-r','linewidth',1.5)
L2_r2r1 = sqrt(sum((r2(:,3)-r1(:,3)).^2));

load test3_cost_r2_qoqa.mat tspan r_tip
r3 = r_tip;
plot(tspan,r3(:,3),'--xb','linewidth',1.5)
L2_r3r1 = sqrt(sum((r3(:,3)-r1(:,3)).^2));

load test3_cost_r2_qo.mat tspan r_tip
r4 = r_tip;
plot(tspan,r4(:,3),'-r','linewidth',1.5)
L2_r4r1 = sqrt(sum((r4(:,3)-r1(:,3)).^2));
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$z_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$ - $q_0,q_a$','FOM: $n=10$ - $q_0$','ROM: $r=2$ - $q_0,q_a$','ROM: $r=2$ - $q_0$','FontSize',14,'Interpreter','latex')



load test3_10sin20t.mat tspan r_tip % FOM simulation data
x_true = csvread('test3_ax_data.csv'); % Ground truth data (BROWN,2003)
z_true = csvread('test3_az_data.csv'); % Ground truth data (BROWN,2003)

figure
plot(tspan,r_tip(:,1),'-r','linewidth',1.5)
hold on
plot(tspan(2:end),x_true+1,'--k','linewidth',1.5)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$x_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 1000, 250])
figure
plot(tspan,r_tip(:,3),'-r','linewidth',1.5)
hold on
plot(tspan(2:end),z_true,'--k','linewidth',1.5)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$z_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 1000, 250])


load test3_10sin50t.mat tspan r_tip
x_true = csvread('test3_bx_data.csv'); % Ground truth data (BROWN,2003)
z_true = csvread('test3_bz_data.csv'); % Ground truth data (BROWN,2003)

figure
plot(tspan(1:10:end),r_tip(1:10:end,1),'-r','linewidth',2)
hold on
plot(0.005:0.005:2,x_true+1,'--k','linewidth',1.5)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$x_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 1000, 250])
figure
plot(tspan(1:10:end),r_tip(1:10:end,3),'-r','linewidth',2)
hold on
plot(0.005:0.005:2,z_true,'--k','linewidth',1.5)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$z_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 1000, 250])


load test3_5sin55-6t.mat tspan r_tip
x_true = csvread('test3_cx_data.csv'); % Ground truth data (BROWN,2003)
z_true = csvread('test3_cz_data.csv'); % Ground truth data (BROWN,2003)

figure
plot(tspan(1:100:end),r_tip(1:100:end,1),'-r','linewidth',2)
hold on
plot(0:0.005:1,x_true+1,'--k','linewidth',1.5)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$x_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 1000, 250])
figure
plot(tspan(1:100:end),r_tip(1:100:end,3),'-r','linewidth',2)
hold on
plot(0:0.005:1,z_true,'--k','linewidth',1.5)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$z_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=10$','Ground Truth','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 1000, 250])


global beam
element_positionCM(Xeq,DoF);
r = [beam(1).r0(1,:);beam(1).r1];
rCM = beam(1).rCM;
r(:,2) = [rCM(1,2);rCM(:,2)];
r_tip = r(end,:);

figure
plot(r(:,1),r(:,3),'-or','linewidth',2)
xlabel('$x_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
ylabel('$z_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
set(gcf, 'Position',  [500, 500, 1000, 250])


%% Highly Flexible Wing - y_ea=0; y_cg=0.267c - Equilibrium:
load testEq_HFW_FOM.mat r_tip ang_tip t n
figure
plot(n,(r_tip(:,1)-r_tip(1,1))/(r_tip(end,1)-r_tip(1,1)),'-b','linewidth',1.5)
hold on
plot(n,(r_tip(:,3)-r_tip(1,3))/(r_tip(end,3)-r_tip(1,3)),'-r','linewidth',1.5)
plot(n,(ang_tip(:)-ang_tip(1))/(ang_tip(end)-ang_tip(1)),'color',[0 0.5 0],'linewidth',1.5)
plot(n,t/t(end),'-k','linewidth',1.5)
xlabel('$n$','FontSize',16,'Interpreter','latex')
ylabel('Normalized Quantities','FontSize',16,'Interpreter','latex')
grid on
legend('$x_{\textrm{tip}}$','$z_{\textrm{tip}}$','$\theta_{\textrm{tip}}$','t','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 1200, 450])

load testEq_HFW_ROM.mat r_tip ang_tip t r
r = 1:28;
figure
plot(r,(r_tip(:,1)-r_tip(1,1))/(r_tip(end,1)-r_tip(1,1)),'-b','linewidth',1.5)
hold on
plot(r,(r_tip(:,3)-r_tip(1,3))/(r_tip(end,3)-r_tip(1,3)),'-r','linewidth',1.5)
plot(r,(ang_tip(:)-ang_tip(1))/(ang_tip(end)-ang_tip(1)),'color',[0 0.5 0],'linewidth',1.5)
plot(r,t/t(end),'-k','linewidth',1.5)
xlabel('$r$','FontSize',16,'Interpreter','latex')
ylabel('Normalized Quantities','FontSize',16,'Interpreter','latex')
grid on
legend('$x_{\textrm{tip}}$','$z_{\textrm{tip}}$','$\theta_{\textrm{tip}}$','t','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 1200, 450])


for i=2:length(r)
    d_x(i) = r_tip(i,1)-r_tip(i-1,1);
    d_z(i) = r_tip(i,3)-r_tip(i-1,3);
    ang(i) = ang_tip(i)-ang_tip(i-1);
end


%% Highly Flexible Wing - y_ea=0; y_cg=0.267c - Dynamics:
load testDin_HFW_FOM.mat r_tip ang_tip tspan
r_tip_FOM = r_tip;
ang_tip_FOM = ang_tip;
load testDin_HFW_ROM.mat r_tip ang_tip tspan
figure
plot(tspan,r_tip_FOM(:,1),'-r','linewidth',2)
hold on 
plot(tspan,r_tip(:,1),'--b','linewidth',2)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$x_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=28$','ROM: $r=12$','FontSize',14,'Interpreter','latex')
ylim([11 20])
set(gcf, 'Position',  [500, 500, 500, 250])
figure
plot(tspan,r_tip_FOM(:,2),'-r','linewidth',2)
hold on 
plot(tspan,r_tip(:,2),'--b','linewidth',2)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$y_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=28$','ROM: $r=12$','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 500, 250])
figure
plot(tspan,r_tip_FOM(:,3),'-r','linewidth',2)
hold on 
plot(tspan,r_tip(:,3),'--b','linewidth',2)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$z_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=28$','ROM: $r=12$','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 500, 250])
figure
plot(tspan,ang_tip_FOM,'-r','linewidth',2)
hold on 
plot(tspan,ang_tip,'--b','linewidth',2)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$\theta_{\textrm{tip}}$ [deg]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM: $n=28$','ROM: $r=12$','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 500, 250])


figure
x_vec = 0:0.01:16;
plot(x_vec,Fy(5,x_vec),'-r','linewidth',2)
hold on 
plot(x_vec,Fz(5,x_vec),'-b','linewidth',2)
xlabel('Half-span position [m]','FontSize',16,'Interpreter','latex')
ylabel('Force [N]','FontSize',16,'Interpreter','latex')
grid on
legend('$F_y$','$F_z$','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 1000, 250])


%% MLP and DHNN Tests:

load test.mat
q_FOM = X_FOM(:,size(X_FOM,2)/2+1:end);
q_MLP = X_MLP(:,size(X_MLP,2)/2+1:end);
q_DHNN = X_DHNN(:,size(X_DHNN,2)/2+1:end);
p_FOM = X_FOM(:,1:size(X_FOM,2)/2);
p_MLP = X_MLP(:,1:size(X_MLP,2)/2);
p_DHNN = X_DHNN(:,1:size(X_DHNN,2)/2);

global beam
for i = 1:size(q_FOM,1)
    element_positionCM(q_FOM(i,:),DoF);
    r = [beam(1).r0(1,:);beam(1).r1];
    rCM = beam(1).rCM;
    r(:,2) = [rCM(1,2);rCM(:,2)];
    r_tip_FOM(i,:) = r(end,:);
end
for i = 1:size(q_MLP,1)
    element_positionCM(q_MLP(i,:),DoF);
    r = [beam(1).r0(1,:);beam(1).r1];
    rCM = beam(1).rCM;
    r(:,2) = [rCM(1,2);rCM(:,2)];
    r_tip_MLP(i,:) = r(end,:);
end
for i = 1:size(q_DHNN,1)
    element_positionCM(q_DHNN(i,:),DoF);
    r = [beam(1).r0(1,:);beam(1).r1];
    rCM = beam(1).rCM;
    r(:,2) = [rCM(1,2);rCM(:,2)];
    r_tip_DHNN(i,:) = r(end,:);
end

figure
plot(tvec,r_tip_FOM(:,3),'--k','linewidth',2)
hold on 
plot(tvec,r_tip_MLP(:,3),'-r','linewidth',2)
plot(tvec,r_tip_DHNN(:,3),'-b','linewidth',2)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$z_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM','MLP', 'D-HNN','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 700, 550])
figure
plot(tvec,r_tip_FOM(:,1),'--k','linewidth',2)
hold on 
plot(tvec,r_tip_MLP(:,1),'-r','linewidth',2)
plot(tvec,r_tip_DHNN(:,1),'-b','linewidth',2)
xlabel('$t$ [s]','FontSize',16,'Interpreter','latex')
ylabel('$x_{\textrm{tip}}$ [m]','FontSize',16,'Interpreter','latex')
grid on
legend('FOM','MLP', 'D-HNN','FontSize',14,'Interpreter','latex')
set(gcf, 'Position',  [500, 500, 700, 550])