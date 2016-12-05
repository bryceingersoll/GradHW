clear; clc; close all;

x = 0.2;
u_max = 0.35;
rho = 1.205;
nu = 15.11*10^-6;
mu = rho*nu;

% laminar jet
J = u_max*nu*x*8*pi/3;

width = 2*5.269*nu*x*J^(-0.5);

Q = 8*pi*nu*x;

m_dot = Q*rho;

%turbulent jet
