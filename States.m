%%%%%%%%% car?like robot with constant speed %%%
clc
clear all
close all
syms x1 x2 x3 x4 x5 x6 v L 
X=[x1 x2 x3 x4 x5 x6];

syms a0 a1 a2 a3 k1 k2 k3 k5 k6

fid = fopen('states.txt', 'wt');
%% System
fx=[(v+x5)*cos(x3);(v+x5)*sin(x3);((v+x5)/L)*tan(x4);0;x6;0];
g1x=[0;0;0;0;0;1];
g2x = [0;0;0;1;0;0];
%% Path
%S=x2?cos(x1); %% Sinusoidal path
S=(x1)^2+(x2)^2-(1^2);
fprintf(fid, '%s\n','xi_1=', char(S));
%% First Lie Derivative
pd=jacobian(S,X);
LfS=pd*fx;
Lg1S=pd*g1x;
Lg2S=pd*g2x;
S_dot=LfS+Lg1S+Lg2S;
S_dot_simplify=simplify(S_dot);
fprintf(fid, '%s\n','xi_2=', char(S_dot));
pretty(S_dot_simplify);
%% Second Lie derivative
pd=jacobian(S_dot,X);
Lf2S=pd*fx;
fprintf(fid, '%s\n','Lf2S=', char(Lf2S));

Lg1LfS=pd*g1x;
Lg2LfS=pd*g2x;
S_ddot=simplify(Lf2S+Lg1LfS+Lg2LfS);
fprintf(fid, '%s\n','xi_3=', char(S_ddot));
% % %
%% Third Lie Derivative

pd=[diff(S_ddot,x1) diff(S_ddot,x2) diff(S_ddot,x3) diff(S_ddot,x4) diff(S_ddot,x5) diff(S_ddot,x6)];
Lf3S=pd*fx;
fprintf(fid, '%s\n', 'Lf3S = ', char(Lf3S));
Lg1Lf2S=pd*g1x;
fprintf(fid, '%s\n', 'Lg1Lf2S = ', char(Lg1Lf2S));
Lg2Lf2S=pd*g2x;
fprintf(fid, '%s\n', 'Lg2Lf2S = ', char(Lg2Lf2S));
S_d3dot=simplify(Lf3S);

%% Eta states
P=atan(x2/x1);
fprintf(fid, '%s\n', 'eta_1 = ', char(P));
%% First Lie Derivative
pd=jacobian(P,X);
LfP=pd*fx;
Lg1P=pd*g1x;
Lg2P=pd*g2x;
P_dot=simplify(LfP+Lg1P+Lg2P);
eta2=P_dot;
fprintf(fid, '%s\n','eta_2=',char(eta2));
%% Second Lie derivative
pd=jacobian(P_dot,X);
Lf2P=pd*fx;
Lg1LfP=pd*g1x;
Lg2LfP=pd*g2x;
P_ddot=simplify(Lf2P+Lg1LfP+Lg2LfP);
eta3=P_ddot;
fprintf(fid, '%s\n','eta_3 = ', char(eta3));
%% Third Derivative
pd=jacobian(P_ddot,X);
Lf3P=pd*fx;
fprintf(fid, '%s\n', 'Lf3P = ' ,char(Lf3P));
Lg1Lf2P=simplify(pd*g1x);
fprintf(fid, '%s\n', 'Lg1Lf2P = ' ,char(Lg1Lf2P));
Lg2Lf2P=simplify(pd*g2x);
fprintf(fid, '%s\n', 'Lg2Lf2P = ' ,char(Lg2Lf2P));
P_d3dot=simplify(Lf3S);

d21=Lg1Lf2P;
fprintf(fid, '%s\n', 'd21 = ' ,char(d21));
d22=Lg2Lf2P;
fprintf(fid, '%s\n', 'd22 = ' ,char(d22));
d11=Lg1Lf2S;
fprintf(fid, '%s\n', 'd11 = ' ,char(d11));
d12=Lg2Lf2S;
fprintf(fid, '%s\n', 'd12 = ' ,char(d12));


D=[d21 d22;d11 d12];

M=simplify(inv(D))
fprintf(fid, '%s\n', 'M = ' ,char(M));

fclose(fid);
