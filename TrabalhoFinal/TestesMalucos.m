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
 ell = 2;
 I = (m*ell^2)/3;
 b = 0.005;
 g = 9.8;

 constanteA = I + m*ell^2; 
 ml = m * ell; %4
 m2l2 = m^2*ell^2; %16
 Mm = m + M; %12
 
 %%fun��es A
 f1 = xP;
 f2 = ((u - b*xP)/Mm - (m2l2*g*theta)/(Mm*constanteA))/(1-(m2l2)/(Mm*constanteA));
 f3 = thetaP;
 f4 = ml*(g*theta-(u+b*xP)/Mm)/(constanteA + m2l2/Mm);

Aaux = [ f1; f2; f3; f4 ];

Jaux = jacobian(Aaux, [x, xP, theta, thetaP ]);

B = [ diff(f1, u); diff(f2, u);diff(f3, u); diff(f4, u) ]; 
C = [ 1 0 0 0; 0 0 1 0 ];
rank(C)
D = [ 0; 0];

Ctrb = ctrb(Jaux, B);
Obsv = obsv(Jaux, C);

rank(Ctrb);
rank(Obsv);


estados = {'x' 'xP' 'theta' 'thetaP'};
entradas = {'u'};
saidas = {'x'; 'theta'};

espacoDeEstados = ss(double(Jaux),double(B),double(C),double(D),'statename',estados,'inputname',entradas,'outputname',saidas);

H = tf(espacoDeEstados);

F = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; -1296 -864 -216 -24];
Lhat = [ 1 0; 0 1; 0 0; 0 0];
rank(ctrb(F,Lhat))
T = lyap(-double(F), double(Jaux), -double(Lhat * C));
L = inv(T)*Lhat;