%%parametros do sistema
syms theta;
syms thetaP;
syms theta2P;
syms x;
syms xP;
syms x2P;
syms u;

%%variaveis do processo
syms m;
syms M;
syms l;
syms I;
syms b;
syms g;

Aaux = [ thetaP; (m*g*l*theta - m*l*x2P)/(I+m*l^2); xP; (u - b*xP)/(M+m) ]

Jaux = jacobian(Aaux, [theta, theta2P, xP, x2P])

rank(jacobian(Aaux))

ctrb(jacobian(Aaux))