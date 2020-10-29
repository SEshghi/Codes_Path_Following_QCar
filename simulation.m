%Path invarient Control - Mltiple Inpts
clc
clear
close all
%% =========== Set the paramters =======
dt=0.01; % Sampling Time
tfinal=100; % final simlation time
time = 0:dt:tfinal;

L=0.3;
d1=0.1;
%% =========== Initialize the states =======
x1 = 1.5; % initilize the position in x direction       - represented as x_1
x2 = 1.2; % initilize the position in y direction       - represented as x_2
x3 = pi/2; % initilize the heading                       - represented as x_3
x4 = 0.1; % initilize the derivative of linear velocity - represented as x_4
x5 = 0.1; 
x6 = 0.3;

k1 = 200;
k2 = 100;
k3 = 60;

k5 = 50;
k6 = 30;
%% =========== velocity profile =======
v = 0.5;   % constant linear velocity                   
%=====================================

x_1_plot = zeros(1,length(time)-1);
x_2_plot = zeros(1,length(time)-1);
x_3_plot = zeros(1,length(time)-1);
x_4_plot = zeros(1,length(time)-1);
x_5_plot = zeros(1,length(time)-1);
x_6_plot = zeros(1,length(time)-1);
eta_1_plot = zeros(1,length(time)-1);
eta_2_plot = zeros(1,length(time)-1);
eta_3_plot = zeros(1,length(time)-1);
xi_1_plot = zeros(1,length(time)-1);
xi_2_plot = zeros(1,length(time)-1);
xi_3_plot = zeros(1,length(time)-1);
xi_2_plot = zeros(1,length(time)-1);
%% =========== The main loop ==========
for k = 1:tfinal/dt

% % Desired velocity profile
%     if k>=0 && k <10
%         eta_ref_2 = -0.5;
%     
%     else 
%         eta_ref_2 = -1.5;
%     end


%   xi states
xi_1 = (x1)^2 + (x2)^2 - 1;
xi_2 = 2*x1*cos(x3)*(v + x5) + 2*x2*sin(x3)*(v + x5);
xi_3 = 2*cos(x3)^2*(v + x5)^2 + 2*sin(x3)^2*(v + x5)^2 + x6*(2*x1*cos(x3) + 2*x2*sin(x3)) + (tan(x4)*(2*v + 2*x5)*(v + x5)*(x2*cos(x3) - x1*sin(x3)))/L;


  
Lf3S = x6*(2*cos(x3)^2*(2*v + 2*x5) + 2*sin(x3)^2*(2*v + 2*x5) + (2*tan(x4)*(v + x5)*(x2*cos(x3) - x1*sin(x3)))/L + (tan(x4)*(2*v + 2*x5)*(x2*cos(x3) - x1*sin(x3)))/L) + cos(x3)*(v + x5)*(2*x6*cos(x3) - (sin(x3)*tan(x4)*(2*v + 2*x5)*(v + x5))/L) + sin(x3)*(v + x5)*(2*x6*sin(x3) + (cos(x3)*tan(x4)*(2*v + 2*x5)*(v + x5))/L) + (tan(x4)*(x6*(2*x2*cos(x3) - 2*x1*sin(x3)) - (tan(x4)*(2*v + 2*x5)*(v + x5)*(x1*cos(x3) + x2*sin(x3)))/L)*(v + x5))/L;

     
   
% eta states
eta_1= atan(x2/x1);
eta_2= -((v + x5)*(x2*cos(x3) - x1*sin(x3)))/(x1^2 + x2^2);
eta_3 = (cos(x3)*(v + x5)^2*(x2^2*sin(x3) - x1^2*sin(x3) + 2*x1*x2*cos(x3)))/(x1^2 + x2^2)^2 - (x6*(x2*cos(x3) - x1*sin(x3)))/(x1^2 + x2^2) - (sin(x3)*(v + x5)^2*(x1^2*cos(x3) - x2^2*cos(x3) + 2*x1*x2*sin(x3)))/(x1^2 + x2^2)^2 + (tan(x4)*(v + x5)^2*(x1*cos(x3) + x2*sin(x3)))/(L*(x1^2 + x2^2));
  

Lf3P = x6*((cos(x3)*(2*v + 2*x5)*(x2^2*sin(x3) - x1^2*sin(x3) + 2*x1*x2*cos(x3)))/(x1^2 + x2^2)^2 - (sin(x3)*(2*v + 2*x5)*(x1^2*cos(x3) - x2^2*cos(x3) + 2*x1*x2*sin(x3)))/(x1^2 + x2^2)^2 + (tan(x4)*(2*v + 2*x5)*(x1*cos(x3) + x2*sin(x3)))/(L*(x1^2 + x2^2))) + cos(x3)*(v + x5)*((x6*sin(x3))/(x1^2 + x2^2) + (cos(x3)*(v + x5)^2*(2*x2*cos(x3) - 2*x1*sin(x3)))/(x1^2 + x2^2)^2 - (sin(x3)*(v + x5)^2*(2*x1*cos(x3) + 2*x2*sin(x3)))/(x1^2 + x2^2)^2 + (2*x1*x6*(x2*cos(x3) - x1*sin(x3)))/(x1^2 + x2^2)^2 - (4*x1*cos(x3)*(v + x5)^2*(x2^2*sin(x3) - x1^2*sin(x3) + 2*x1*x2*cos(x3)))/(x1^2 + x2^2)^3 + (4*x1*sin(x3)*(v + x5)^2*(x1^2*cos(x3) - x2^2*cos(x3) + 2*x1*x2*sin(x3)))/(x1^2 + x2^2)^3 + (cos(x3)*tan(x4)*(v + x5)^2)/(L*(x1^2 + x2^2)) - (2*x1*tan(x4)*(v + x5)^2*(x1*cos(x3) + x2*sin(x3)))/(L*(x1^2 + x2^2)^2)) + sin(x3)*(v + x5)*((cos(x3)*(v + x5)^2*(2*x1*cos(x3) + 2*x2*sin(x3)))/(x1^2 + x2^2)^2 - (x6*cos(x3))/(x1^2 + x2^2) + (sin(x3)*(v + x5)^2*(2*x2*cos(x3) - 2*x1*sin(x3)))/(x1^2 + x2^2)^2 + (2*x2*x6*(x2*cos(x3) - x1*sin(x3)))/(x1^2 + x2^2)^2 - (4*x2*cos(x3)*(v + x5)^2*(x2^2*sin(x3) - x1^2*sin(x3) + 2*x1*x2*cos(x3)))/(x1^2 + x2^2)^3 + (4*x2*sin(x3)*(v + x5)^2*(x1^2*cos(x3) - x2^2*cos(x3) + 2*x1*x2*sin(x3)))/(x1^2 + x2^2)^3 + (sin(x3)*tan(x4)*(v + x5)^2)/(L*(x1^2 + x2^2)) - (2*x2*tan(x4)*(v + x5)^2*(x1*cos(x3) + x2*sin(x3)))/(L*(x1^2 + x2^2)^2)) + (tan(x4)*(v + x5)*((x6*(x1*cos(x3) + x2*sin(x3)))/(x1^2 + x2^2) - (2*cos(x3)*(v + x5)^2*(x1^2*cos(x3) - x2^2*cos(x3) + 2*x1*x2*sin(x3)))/(x1^2 + x2^2)^2 - (2*sin(x3)*(v + x5)^2*(x2^2*sin(x3) - x1^2*sin(x3) + 2*x1*x2*cos(x3)))/(x1^2 + x2^2)^2 + (tan(x4)*(v + x5)^2*(x2*cos(x3) - x1*sin(x3)))/(L*(x1^2 + x2^2))))/L;



Lg1Lf2P = -(x2*cos(x3) - x1*sin(x3))/(x1^2 + x2^2);
Lg2Lf2P = ((v + x5)^2*(x1*cos(x3) + x2*sin(x3)))/(L*cos(x4)^2*(x1^2 + x2^2));

Lg1Lf2S = 2*x1*cos(x3) + 2*x2*sin(x3);
Lg2Lf2S = ((2*v + 2*x5)*(v + x5)*(tan(x4)^2 + 1)*(x2*cos(x3) - x1*sin(x3)))/L;


D=[Lg1Lf2P Lg2Lf2P;Lg1Lf2S Lg2Lf2S];

M=inv(D);

v1_fb1 = -Lf3P -k5*(eta_2-0.5)-k6*(eta_3-0);
v2_fb1 = -Lf3S -k1*xi_1-k2*xi_2-k3*xi_3;   

   U=M*[v1_fb1;v2_fb1];
   
   u1=U(1);
   u2=U(2);
   

  
% inpt commands
vel = x5 + v;
% Kinematics with dynamic extension
    x1 = vel*cos(x3)*dt+x1; % calclating x
    x2 = vel*sin(x3)*dt+x2; % calclating y
    x3 = (vel/L)*tan(x4)*dt+x3;
    x4 = u2*dt+x4; % calclating theta // _2 = w anglar velocity
    x6 = u1*dt+x6;
    x5 = x6*dt + x5;
    
    
    x_1_plot(k) = x1;
    x_2_plot(k) = x2;
    x_3_plot(k) = x3;
    x_4_plot(k) = x4;
     x_5_plot(k) = x4;
      x_6_plot(k) = x4;
    eta_1_plot(k) = eta_1;
    eta_2_plot(k) = eta_2;
     eta_2_plot(k) = eta_3;
    xi_1_plot(k) = xi_1;
    xi_2_plot(k) = xi_2;
    xi_3_plot(k) = xi_3;
    u_2_plot(k) = 2;
  
end

%% =========== Plot the reslts =======

figure(1)
plot(x_1_plot,x_2_plot,'b--','LineWidth',2)
xlabel('$x(m)$','FontSize',16,'Interpreter','latex')
ylabel('$y (m)$','FontSize',16,'Interpreter','latex')
% xlim([-1.1 1.1])
% ylim([-1.1 1.1])
% hold on
% h = desired_path(0,0,1);



figure(2)
plot(time(1:end-1),xi_1_plot,'b','LineWidth',2)
xlabel('$t(s)$','FontSize',16,'Interpreter','latex')
ylabel('$\xi_1$','FontSize',16,'Interpreter','latex')

figure(3)
plot(time(1:end-1),xi_2_plot,'b','LineWidth',2)
xlabel('$t(s)$','FontSize',16,'Interpreter','latex')
ylabel('$\xi_2$','FontSize',16,'Interpreter','latex')

figure(4)
plot(time(1:end-1),eta_1_plot,'b','LineWidth',1)
xlabel('$t(s)$','FontSize',16,'Interpreter','latex')
ylabel('$\eta_1$','FontSize',16,'Interpreter','latex')

figure(5)
plot(time(1:end-1),eta_2_plot,'b','LineWidth',1)
xlabel('$t(s)$','FontSize',16,'Interpreter','latex')
ylabel('$\eta_2$','FontSize',16,'Interpreter','latex')


figure(6)
plot(time(1:end-1),eta_2_plot,'b','LineWidth',1)
xlabel('$t(s)$','FontSize',16,'Interpreter','latex')
ylabel('$\eta_2$','FontSize',16,'Interpreter','latex')

figure(7)
plot(time(1:end-1),u_2_plot,'b','LineWidth',1)
xlabel('$t(s)$','FontSize',16,'Interpreter','latex')
ylabel('$\omega$','FontSize',16,'Interpreter','latex')
% %=====================================
