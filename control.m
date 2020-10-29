 %%%%%%%%% car?like robot attached with one trailer system %%%
clc
clear all
close all
syms x1 x2 x3 x4 x5 x6 x7 v L u d1
syms a0 a1 a2 a3 k1 k2
X=[x1 x2 x3 x4 x5 x6 x7];

fid = fopen('control.txt','wt');

%% System
fx=[(v+x6)*cos(x3);(v+x6)*sin(x3);((v+x6)/L)*tan(x4);0;((v+x6)/d1)*sin(x3-x5);x7;0];
g1x=[0;0;0;1;0;0;0];
g2x=[0;0;0;0;0;0;1];
%% Path: Circular Path
S=(x1-L*cos(x5))^2+(x2-L*sin(x5))^2-1;
fprintf(fid, '%s\n','xi_1=', char(S));
%% First Lie Derivative
pd=[diff(S,x1) diff(S,x2) diff(S,x3) diff(S,x4) diff(S,x5) diff(S,x6) ...
diff(S,x7)];
LfS=pd*fx;
Lg1S=pd*g1x*u;
Lg2S=pd*g2x*u;
S_dot=LfS+Lg1S+Lg2S;
S_dot_simplify=simplify(S_dot);
fprintf(fid, '%s\n','xi_2=', char(S_dot));
pretty(S_dot_simplify);
%% Second Lie derivative
pd=[diff(S_dot,x1) diff(S_dot,x2) diff(S_dot,x3) diff(S_dot,x4) ...
diff(S_dot,x5) diff(S_dot,x6) diff(S_dot,x7)];
Lf2S=pd*fx;
fprintf(fid, '%s\n','Lf2S=', char(Lf2S));

Lg1LfS=pd*g1x*u;
Lg2LfS=pd*g2x*u;
S_ddot=simplify(Lf2S+Lg1LfS+Lg2LfS);
fprintf(fid, '%s\n','xi_3=', char(S_ddot));
% % %
%% Third Lie Derivative

pd=[diff(S_ddot,x1) diff(S_ddot,x2) diff(S_ddot,x3) diff(S_ddot,x4) ...
diff(S_ddot,x5) diff(S_ddot,x6) diff(S_ddot,x7)];
Lf3S=pd*fx;
Lg1Lf2S=pd*g1x;
Lg2Lf2S=pd*g2x;
S_d3dot=simplify(Lf3S);

%% Forth Derivative
pd=[diff(S_d3dot,x1) diff(S_d3dot,x2) diff(S_d3dot,x3) ...
diff(S_d3dot,x4) diff(S_d3dot,x5) diff(S_d3dot,x6) diff(S_d3dot,x7)];
Lf4S=pd*fx;
Lg1Lf3S=pd*g1x;
fprintf(fid, '%s\n', 'Lg1Lf3S = ',char(Lg1Lf3S));

Lg2Lf3S=pd*g2x;
fprintf(fid, '%s\n', 'Lg2Lf3S = ' ,char(Lg2Lf3S));

%% Eta states
P=atan( (x2-L*sin(x5)) / (x1-L*cos(x5)) );
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
Lg1Lf2P=simplify(pd*g1x);
fprintf(fid, '%s\n', 'Lg1Lf2P = ' ,char(Lg1Lf2P));

Lg2Lf2P=simplify(pd*g2x);
fprintf(fid, '%s\n', 'Lg2Lf2P = ' ,char(Lg2Lf2P));

P_d3dot=simplify(Lf3S);

%% Fourth Derivative Derivative
pd=jacobian(P_d3dot,X);
Lf4P=pd*fx;
fprintf(fid, '%s\n', 'Lf4P = ' ,char(Lf4P));
Lg1Lf3P=simplify(pd*g1x);
Lg2Lf3P=simplify(pd*g2x);

d21=Lg1Lf2P;
d22=Lg2Lf2P;
d11=Lg1Lf3S;
d12=Lg2Lf3S;


D=[d21 d22;d11 d12];
M=inv(D);
fclose(fid);