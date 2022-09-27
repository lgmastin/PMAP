function [up_t] = get_terminal_velocity(dp,rho_p,rho_f,nu,g,vargin)
% iterate to find the terminal velocity of particle classes

if nargin==6
    up_t = vargin; % initialize with guess velocity
else
    up_t = ones(1,numel(dp)); % otherwise initialize with ones
end

resid = Inf;
tol = 0.01;
iter = 0;
while any(resid>tol)
    up_t_old = up_t;
    Re = up_t.*dp/nu;
    Cd = (Re<=1000).*24./Re.*(1+0.15*Re.^0.687)+(Re>1000)*0.44; % Schiller-Naumann correlation
    up_t = sqrt(4/3*dp./Cd.*rho_p/rho_f*g);
    resid = abs(up_t-up_t_old)./up_t_old;
    iter = iter+1;
end

if any(Re>2e5)
    disp('Warning: Drag catastrophe exceeded')
end
    
end