@@ -0,0 +1,87 @@
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
D = [ 0; 0];

Ctrb = ctrb(Jaux, B);
Obsv = obsv(Jaux, C);

estados = {'x' 'xP' 'theta' 'thetaP'};
entradas = {'u'};
saidas = {'x'; 'theta'};
espacoDeEstados = ss(double(Jaux),double(B),double(C),double(D),'statename',estados,'inputname',entradas,'outputname',saidas);

H = tf(espacoDeEstados);

F = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; -1296 -864 -216 -24];
Lhat = [ 1 0; 0 1; 0 0; 0 0];
T = lyap(-double(F), double(Jaux), -double(Lhat * C));
L = inv(T)*Lhat;

R = [0 0 0 1; 0 1 0 0];
X = [C;R];
Y = inv(X);

Y1 = Y(:,1:2);
Y2 = Y(:,3:4);

Ab = (X*double(Jaux))/X;
Bb = X*double(B);
Cb = double(C)/X;
Db = double(D);

% Capturando as submatrizes de interesse para o pendulo invertido analisado
A11 = Ab(1:2,1:2);
A12 = Ab(1:2,3:4);
A21 = Ab(3:4,1:2);
A22 = Ab(3:4,3:4);
B1 = Bb(1:2,:);
B2 = Bb(3:4,:);
C1 = Cb(:,1:2);
C2 = Cb(:,3:4);

novoLhat = [1 0; 0 0];

T = lyap(-double(F), double(Ab), -double(novoLhat * C1));
L = inv(T)*novoLhat;
Ab - Bb
P = care(Ab, Bb, Cb'*Cb, 1);
KLqr = Bb'*P;
MLqr = KLqr'/(KLqr*KLqr');

syms s;
Cb/(s*eye(4)-(A-Bb*2))*Bb*3