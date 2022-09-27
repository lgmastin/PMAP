function [x_next,y_next,h_next] = advance_RKCK_adapt(x,y,h,tol)
yscale = y;
yscale(5:7) = 1;
yscale(8:end) = max(y(8:end));
S = 0.9; % safety factor for estimating new time step size

while true
    [y_next,yerr] = advance_RKCK(x,y,h);
    if isnan(y_next) % if advancement produces imaginary result, step down h
        h = h/2;
        continue
    end
    max_error=max(abs(yerr./yscale)); % find maximum error, scaled to scale factors
    max_error=max_error/tol;
    if max_error>1
        alpha=max_error^0.25; % adjustment following Dormund, J.R.
        h=S*h/min(alpha,10); % "Numerical Methods for Differential Equations", p. 85
        x_next=x+h;
        continue
    else
        % increase h, but not by more than 5x
        h_next=min(S*h*max_error^-0.2,5*h);
        x_next=x+h;
        break
    end
end

end

%%
function [y_next,y_err] = advance_RKCK(x,y,h)

A2 = 0.2; A3 = 0.3; A4 = 0.6; A5 = 1.0; A6 = 0.875;
B21 = 0.2; B31 = 3/40; B32 = 9/40; B41 = 0.3; B42 = -0.9; B43 = 1.2; B51 = -11/54;
B52 = 2.5; B53 = -70/27; B54 = 35/27; B61 = 1631/55296; B62 = 175/512;
B63 = 575/13824; B64 = 44275/110592; B65 = 253/4096;
C1 = 37/378; C3 = 250/621; C4 = 125/594; C6 = 512/1771;
DC1 = C1-2825/27648; DC3 = C3-18575/48384; DC4 = C4-13525/55296;
DC5 = -277/14336; DC6 = C6-0.25;

y_next = NaN;
y_err = NaN;

dydx = derivs(x,y);
ak1 = h*dydx;
ytmp = y+B21*ak1;
if ~isreal(ytmp)
    return
end
dydx = derivs(x+A2*h,ytmp);
ak2 = h*dydx;
ytmp = y+B31*ak1+B32*ak2;
if ~isreal(ytmp)
    return
end
dydx = derivs(x+A3*h,ytmp);
ak3 = h*dydx;
ytmp = y+B41*ak1+B42*ak2+B43*ak3;
if ~isreal(ytmp)
    return
end
dydx = derivs(x+A4*h,ytmp);
ak4 = h*dydx;
ytmp = y+B51*ak1+B52*ak2+B53*ak3+B54*ak4;
if ~isreal(ytmp)
    return
end
dydx = derivs(x+A5*h,ytmp);
ak5 = h*dydx;
ytmp = y+B61*ak1+B62*ak2+B63*ak3+B64*ak4+B65*ak5;
if ~isreal(ytmp)
    return
end
dydx = derivs(x+A6*h,ytmp);
ak6 = h*dydx;
ytmp = y+C1*ak1+C3*ak3+C4*ak4+C6*ak6;
if ~isreal(ytmp)
    return
end

y_next = ytmp;
y_err = h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6);

end
