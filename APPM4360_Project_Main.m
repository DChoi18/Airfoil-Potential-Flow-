%% APPM 4360 Complex Variables Project
% This project examines potential flow over a Joukowski and Karman-Trefftz
% airfoil
clc; clear; close all
% Author: Derrick Choi
% Collaborators: Aneesh Balla, Zach Berriman-Rozen
%% Complex domain
xmin = -10; xmax = 10;
ymin = -10; ymax = 10;
numx = 200; numy = 200;

xvalues = linspace(xmin,xmax,numx);
yvalues = linspace(ymin,ymax,numy);
complexgrid = complex((ones(length(xvalues),1)*xvalues),(ones(length(yvalues),1)*yvalues)');
%% Cylinder and Airfoil flow
Uinf = 50;
%Circle parameters
center = -0.5+0.1i;
b =2;
R = sqrt((real(center)-real(b))^2+(imag(center)-imag(b))^2);
theta = linspace(0,2*pi,100);
circle = R*exp(theta*1i)+center;
% flow parameters
alpha = deg2rad(15);
circulation = -4*pi*R*Uinf*sin(alpha);
Cylinderflow = @(z) Uinf*(R^2./z)+circulation/(2*pi*1i)*log(z/R);
Jgrid = complexgrid; Jgrid(abs(complexgrid-center)<R-eps) =NaN; %remove points in circle
z = (complexgrid-center)*exp(-1i*alpha); %shifted and rotation of circle
w = KarmanTrefftzTransform(Jgrid,b,1.7); %grid in joukowski transform
airfoil = KarmanTrefftzTransform(circle,b,1.7);
AirfoilFlow = Cylinderflow(z)+Uinf*(complexgrid-center)*exp(-1i*alpha);
%% Plot
figure('Name','Potential flow over Airfoil','Position',[228 61 1150 668])
colormap('default')
subplot(2,2,1)
contour(real((z+center)*exp(1i*alpha)),imag((z+center)*exp(1i*alpha)),real(AirfoilFlow),'LevelList',linspace(min(real(AirfoilFlow),[],'all'),max(real(AirfoilFlow),[],'all'),200))
hold on
plot(real(circle*exp(1i*alpha)),imag(circle*exp(1i*alpha)),'k','LineWidth',1.5)
title('\textbf{Equipotential lines}','Interpreter','latex','FontSize',12)
xlabel('Re(z)','Interpreter','latex','FontSize',12)
ylabel('Im(z)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',12)

subplot(2,2,2)
contour(real((z+center)*exp(1i*alpha)),imag((z+center)*exp(1i*alpha)),imag(AirfoilFlow),'LevelList',linspace(min(imag(AirfoilFlow),[],'all'),max(imag(AirfoilFlow),[],'all'),200));
hold on
plot(real(circle*exp(1i*alpha)),imag(circle*exp(1i*alpha)),'k','LineWidth',1.5)
title('\textbf{Streamlines}','Interpreter','latex','FontSize',12)
xlabel('Re(z)','Interpreter','latex','FontSize',12)
ylabel('Im(z)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',12)

subplot(2,2,3)
contour(real(w),imag(w),real(AirfoilFlow),'LevelList',linspace(min(real(AirfoilFlow),[],'all'),max(real(AirfoilFlow),[],'all'),200));
hold on
plot(real(airfoil),imag(airfoil),'k','LineWidth',1.5)
title('\textbf{Equipotential lines}','Interpreter','latex','FontSize',12)
xlabel('Re(w)','Interpreter','latex','FontSize',12)
ylabel('Im(w)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',12)

subplot(2,2,4)
contour(real(w),imag(w),imag(AirfoilFlow),'LevelList',linspace(min(imag(AirfoilFlow),[],'all'),max(imag(AirfoilFlow),[],'all'),200));
hold on
plot(real(airfoil),imag(airfoil),'k','LineWidth',1.5)
title('\textbf{Streamlines}','Interpreter','latex','FontSize',12)
xlabel('Re(w)','Interpreter','latex','FontSize',12)
ylabel('Im(w)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',12)
sgtitle('\textbf{Flow over K\''arm\''an-Trefftz Airfoil}','Interpreter','latex','FontSize',16)
%% Playing with Joukowski Tranformations
clear;clc;close all
center = -.5+0.25i;
b = 2;
theta = linspace(0,2*pi);
r = sqrt((real(center)-real(b))^2+(imag(center)-imag(b))^2);
z = r*exp(1i*theta)+center;
w = JoukowskiTransform(z,b,'forward');
%Plot
figure('Name','Joukowski Transformation Cambered Airfoil','Position',[400 100 800 600])
p1 = plot(real(z),imag(z),'LineWidth',1.5);
hold on
plot(real(center),imag(center),'.k','MarkerSize',10)
p2 = plot(real(w),imag(w),'r','LineWidth',1.5);
xlabel('Real axis','Interpreter','latex','FontSize',14);
ylabel('Imaginary axis','Interpreter','latex','FontSize',14);
set(gca,'TickLabelInterpreter','latex','FontSize',12,'XAxisLocation','origin','YAxisLocation','origin')
legend([p1,p2],'Curve in $z$-plane','Curve in $w$-plane','Interpreter','latex','Fontsize',12,'Location','northeast')
title('\textbf{Joukowski Transformation to Cambered Airfoil}','Interpreter','latex','FontSize',16)
grid on
ylim([-2*b 2*b])
set(gcf,'Color','w','PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r600',['./images/' strrep(get(gcf,'Name'),' ','') '.png'])

%% Playing with Karman Trefftz Transforms
clear;clc;close all
center = -1+.1i;
b = 2;
theta = linspace(0,2*pi,200);
r = sqrt((real(center)-real(b))^2+(imag(center)-imag(b))^2);
z = r*exp(1i*theta)+center;
w = KarmanTrefftzTransform(z,b,2.5);
%Plot the transformation
figure('Name','Karman Trefftz Transformation k2_5','Position',[400 100 800 600])
p1 = plot(real(z),imag(z),'LineWidth',1.5);
hold on
plot(real(center),imag(center),'.k','MarkerSize',10)
p2 = plot(real(w),imag(w),'r','LineWidth',1.5);
xlabel('Real axis','Interpreter','latex','FontSize',14);
ylabel('Imaginary axis','Interpreter','latex','FontSize',14);
set(gca,'TickLabelInterpreter','latex','FontSize',12,'XAxisLocation','origin','YAxisLocation','origin')
legend([p1,p2],'Curve in $z$-plane','Curve in $w$-plane','Interpreter','latex','Fontsize',12,'Location','northeast')
title('\textbf{K\''arm\''an-Trefftz Transformation k = 2.5}','Interpreter','latex','FontSize',16)
ylim([-2*b 2*b])
grid on
set(gcf,'Color','w','PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r600',['./images/' strrep(get(gcf,'Name'),' ','') '.png'])
%% Foward implementation of flow over Joukowski airfoil
xmin = -10; xmax = 10;
ymin = -10; ymax = 10;
numx = 200; numy = 200;

xvalues = linspace(xmin,xmax,numx);
yvalues = linspace(ymin,ymax,numy);

z = complex((ones(length(xvalues),1)*xvalues),(ones(length(yvalues),1)*yvalues)');
Uinf = 50;
center = -0.1+0.1i;
b =2;
alpha = deg2rad(10);
R = sqrt((real(center)-real(b))^2+(imag(center)-imag(b))^2);
z(abs(z-center)<R) = NaN;
z = JoukowskiTransform(z,b,'forward');
circle = R*exp(1i*linspace(0,2*pi))+center;
airfoil = JoukowskiTransform(circle,b,'forward');
Gamma= 4*pi*R*Uinf*b*sin(alpha);
AirfoilFlow = @(z) Uinf*exp(-1i*alpha)*(1/2*(z+sqrt(z.^2-4*b^2))-center)+...
                   Uinf*R^2./(1/2*(z+sqrt(z.^2-4*b^2))-center)+Gamma/(2*pi*1i)*...
                   log((1/2*(z+sqrt(z.^2-4*b^2))-center)/R);
PHI = real(AirfoilFlow(z));
PSI = imag(AirfoilFlow(z));
figure('Position',[400 100 1200 600])
subplot(1,2,1)
contour(real(z),imag(z),PHI,linspace(min(PHI,[],'all')/2,max(PHI,[],'all')/2,50))
hold on
plot(real(airfoil),imag(airfoil))
xlabel('Re(z)')
ylabel('Im(z)')
subplot(1,2,2)
contour(real(z),imag(z),PSI,linspace(min(PSI,[],'all')/2,max(PSI,[],'all')/2,50))
xlabel('Re(z)')
ylabel('Im(z)')
hold on
plot(real(airfoil),imag(airfoil))
%% Simple conformal Maps
clc; clear; close all
%rotation
xmin = -5; xmax = 5;
ymin = -5; ymax = 5;
numx = 50; numy = 50;

xvalues = linspace(xmin,xmax,numx);
yvalues = linspace(ymin,ymax,numy);

z = complex((ones(length(xvalues),1)*xvalues),(ones(length(yvalues),1)*yvalues)');
Uinf = 50;
alpha = deg2rad(10);

UniformFlow = Uinf*z;
Angledflow = exp(-1i*alpha)*z;
figure('Name','Rotation','Position',[250 100 1200 450])
subplot(1,2,1)
contour(real(z),imag(z),real(UniformFlow),linspace(min(real(UniformFlow),[],'all'),max(real(UniformFlow),[],'all'),50),'Color','k')
hold on
contour(real(z),imag(z),imag(UniformFlow),linspace(min(imag(UniformFlow),[],'all'),max(imag(UniformFlow),[],'all'),50),'Color','k')
title('\textbf{z-plane}','Interpreter','latex','FontSize',12)
xlabel('Re(z)','Interpreter','latex','FontSize',12)
ylabel('Im(z)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',12)

subplot(1,2,2)
contour(real(z),imag(z),real(Angledflow),linspace(min(real(Angledflow),[],'all'),max(real(Angledflow),[],'all'),50),'Color','k')
hold on
contour(real(z),imag(z),imag(Angledflow),linspace(min(real(Angledflow),[],'all'),max(real(Angledflow),[],'all'),50),'Color','k')
title('\textbf{w-plane}','Interpreter','latex','FontSize',12)
xlabel('Re(w)','Interpreter','latex','FontSize',12)
ylabel('Im(w)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',12)
sgtitle('\textbf{Pure Rotation}','Interpreter','latex','FontSize',16)

%Translation
circle = 2*exp(1i*linspace(0,2*pi));
move_circle = circle+0.5+1i;

figure('Name','Translation','Position',[400 100 800 600])
p1 = plot(real(circle),imag(circle),'LineWidth',1.5);
hold on
p2 = plot(real(move_circle),imag(move_circle),'LineWidth',1.5,'Color','r');
title('\textbf{Translation}','Interpreter','latex','FontSize',12)
xlabel('Real axis','Interpreter','latex','FontSize',12)
ylabel('Imaginary axis','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',12,'XAxisLocation','origin','YAxisLocation','origin')
legend([p1,p2],'Curve in $z$-plane','Curve in $w$-plane','Interpreter','latex','Fontsize',12,'Location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',12)
grid on

%Dilation
c = 2;
newcirc = c*circle;

figure('Name','dilation','Position',[400 100 800 600])
p1 = plot(real(circle),imag(circle),'LineWidth',1.5);
hold on
p2 = plot(real(newcirc),imag(newcirc),'LineWidth',1.5,'Color','r');
title('\textbf{Dilation}','Interpreter','latex','FontSize',12)
xlabel('Real axis','Interpreter','latex','FontSize',12)
ylabel('Imaginary axis','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',12,'XAxisLocation','origin','YAxisLocation','origin')
legend([p1,p2],'Curve in $z$-plane','Curve in $w$-plane','Interpreter','latex','Fontsize',12,'Location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',12)
grid on

for i = 1:get(gcf,'Number')
    figure(i)
    set(gcf,'Color','w','PaperPositionMode','auto','InvertHardCopy','off')
    print(gcf,'-dpng','-r600',['./images/' strrep(get(gcf,'Name'),' ','') '.png'])
end