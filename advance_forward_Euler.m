function [x_next,y_next,h_next] = advance_forward_Euler(x,y,h)

dydx = derivs(x,y);

h_next = h;
x_next = x+h;
y_next = y+h*dydx;

end
