function h = h_i(T)
% Returns enthalpy of water ice at specified temperature
%
% McBride, B. J. (1993). Coefficients for calculating thermodynamic and 
% transport properties of individual species (Vol. 4513). NASA Langley 
% Research Center.

global R mw_H2O

% check if temperature is within limits
if max(T)>273.15 || min(T)<200
%     disp('Warning: temperature is outside valid range for enthalpy calculation')
end

h = R/mw_H2O*T.*(5.29677970e0-6.75749247e-2*T/2+5.16942109e-4*T.^2/3- ...
    1.43853360e-6*T.^3/4+1.52564794e-9*T.^4/5-3.62266557e4./T);

end