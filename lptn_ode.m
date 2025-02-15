function dT = lptn_ode(t,T)
%LPTN ODE ODE function of LPTN transient
%   Detailed explanation goes here
global A B 
dT = zeros(7,1);

dT = A*T + B;

end
