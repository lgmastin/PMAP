function h = h_v(T)
% Returns enthalpy of water vapor at specified temperature
%
% McBride, B. J. (1993). Coefficients for calculating thermodynamic and 
% transport properties of individual species (Vol. 4513). NASA Langley 
% Research Center.

global R mw_H2O

% check if temperature is within limits
if max(T)>6000 || min(T)<200
    disp('Warning: temperature is outside valid range for enthalpy calculation')
end

h = zeros(1,numel(T));

ind = false(1,numel(T));
ind(T>1000) = true;
Ttmp = T;
T = Ttmp(ind);

h(ind) = R/mw_H2O*T.*(2.67703787e0+2.97318329e-3*T/2-7.73769690e-7*T.^2/3+ ...
    9.44336689e-11*T.^3/4-4.26900959e-15*T.^4/5-2.98858938e4./T);

T = Ttmp(~ind);

h(~ind) = R/mw_H2O*T.*(4.19864056e0-2.03643410e-3*T/2+6.52040211e-6*T.^2/3- ...
    5.48797062e-9*T.^3/4+1.77197817e-12*T.^4/5-3.02937267e4./T);

end