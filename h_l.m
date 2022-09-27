function h = h_l(T)
% Returns enthalpy of liquid water at specified temperature
%
% McBride, B. J. (1993). Coefficients for calculating thermodynamic and 
% transport properties of individual species (Vol. 4513). NASA Langley 
% Research Center.

global R mw_H2O

% check if temperature is within limits
if max(T)>600 || min(T)<273.15
%     disp('Warning: temperature is outside valid range for enthalpy calculation')
end

h = R/mw_H2O*T.*(7.25575005e1-6.62445402e-1*T/2+2.56198746e-3*T.^2/3- ...
    4.36591923e-6*T.^3/4+2.78178981e-9*T.^4/5-4.18865499e4./T);

end