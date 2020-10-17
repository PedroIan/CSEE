close all;
clear all;

%%parametros do sistema
syms theta;
syms thetaP;
syms theta2P;
syms x;
syms xP;
syms x2P;
syms u;

%%variaveis do processo
 m = 2;
 M = 10;
 l = 2;
 I = (m*l^2)/3;
 b = 0.005;
 g = 9.8;

 constanteA = I + m*l^2; 
 ml = m * l; %4
 m2l2 = m^2*l^2; %16
 Mm = m + M; %12
 
 %%fun��es A
 f1 = xP;
 f2 = ((u - b*xP)/Mm - (m2l2*g*theta)/(Mm*constanteA))/(1-(m2l2)/(Mm*constanteA));
 f3 = thetaP;
 f4 = ml*(g*theta-(u+b*xP)/Mm)/(constanteA + m2l2/Mm);

Aaux = [ f1; f2; f3; f4 ];

Jaux = jacobian(Aaux, [theta, thetaP, xP, x2P]);

B = [ diff(f1, u); diff(f2, u);diff(f3, u); diff(f4, u) ]; 
C = [ 1 0 0 0; 0 0 1 0 ]

rank(Jaux)

Ctrb = ctrb(Jaux, B);

PcInv = [ 0 0 0 1; 0.0952 0 1 0; 0 1 0 0; -0.0278 0 0 0 ];
rank(PcInv)
Pc = inv(PcInv);

Ahat = Pc * Jaux * PcInv;
Bhat = Pc * B;
rank(ctrb(Ahat, Bhat))