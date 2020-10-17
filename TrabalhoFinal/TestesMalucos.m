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
 
 %%fun��es A
 f1 = thetaP;
 f2 = (m*g*l*theta - m*l*((u - b*xP)/(M+m)))/(I+m*l^2);
 f3 = xP;
 f4 = (u - b*xP)/(M+m);

Aaux = [ f1; f2; f3; f4 ]

Jaux = jacobian(Aaux, [theta, theta2P, xP, x2P])

B = [ diff(f1, u); diff(f2, u);diff(f3, u); diff(f4, u) ]; 

rank(Jaux)

Ctrb = ctrb(Jaux, B);

PcInv = [ 0 0 0 1; -0.0313 1 0 0; 0 0 1 0; 0.0833 0 0 0 ];
rank(PcInv)
Pc = inv(PcInv);

Ahat = Pc * Jaux * PcInv;
Bhat = Pc * B;
rank(ctrb(Ahat, Bhat))