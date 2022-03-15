airfoil_data = readtable('NACA_0012.txt');
airfoil_data = airfoil_data{:,:};

X3D = repmat(airfoil_data(:,1), 1, 3);
Y3D = repmat(linspace(0, 10, 3), length(airfoil_data(:,1)), 1);
Z3D = repmat(airfoil_data(:,2), 1, 3);
X3D(:,2) = 0.5*X3D(:,2);
X3D(:,3) = 0.25*X3D(:,3);
Z3D(:,2) = 0.5*Z3D(:,2);
Z3D(:,3) = 0.25*Z3D(:,3);
figure
surf(X3D, Y3D, Z3D,'linestyle','none');
hold on
plot3(X3D,Y3D,Z3D,'k')
axis equal;
%set(gca,'Visible','off');

for i = 1:length(airfoil_data(:,1))
pos(i,:) = DCM(3,pi/2)*[airfoil_data(i,1)-0.25;airfoil_data(i,2);0];
end

figure
plot(airfoil_data(:,1)-0.25,airfoil_data(:,2))
hold on
plot(pos(:,1),pos(:,2))
    

