function h = h_a(T)
% Returns enthalpy of dry binary air at specified temperature
%
% McBride, B. J. (1993). Coefficients for calculating thermodynamic and 
% transport properties of individual species (Vol. 4513). NASA Langley 
% Research Center.

global R

% check if temperature is within limits
if max(T)>6000 || min(T)<200
    disp('Warning: temperature is outside valid range for enthalpy calculation')
end

n_O2 = 0.21; % molar fraction of O2
n_N2 = 0.79; % molar fraction of N2

mw_O2 = 32; % molecular weight of O2
mw_N2 = 28; % molecular weight of N2

m_O2 = mw_O2*n_O2/(mw_O2*n_O2+mw_N2*n_N2); % mass fraction of O2 in dry air
m_N2 = mw_N2*n_N2/(mw_O2*n_O2+mw_N2*n_N2); % mass fraction of N2 in dry air

% enthalpies of component species in dry air
h_O2 = zeros(1,numel(T)); h_N2 = h_O2;
ind = false(1,numel(T));
ind(T>1000) = true;
Ttmp = T;
T = Ttmp(ind);

% polynomials in temp range T>1000K
h_O2(ind) = R/mw_O2*T.*(3.66096083e0+6.56365523e-4*T/2-1.41149485e-7*T.^2/3+ ...
    2.06797658e-11*T.^3/4-1.29913248e-15*T.^4/5-1.21597725e3./T);
h_N2(ind) = R/mw_N2*T.*(2.95257626e0+1.39690057e-3*T/2-4.92631691e-7*T.^2/3+ ...
    7.86010367e-11*T.^3/4-4.60755321e-15*T.^4/5-9.23948645e2./T);

T = Ttmp(~ind);

% polynomials in temp range T<=1000K
h_O2(~ind) = R/mw_O2*T.*(3.78245636e0-2.99673416e-3*T/2+9.84730201e-6*T.^2/3- ...  
    9.68129509e-9*T.^3/4+3.24372837e-12*T.^4/5-1.06394356e3./T);
h_N2(~ind) = R/mw_N2*T.*(3.53100528e0-1.23660987e-4*T/2-5.02999437e-7*T.^2/3+ ...
    2.43530612e-9*T.^3/4-1.40881235e-12*T.^4/5-1.04697628e3./T);

T = Ttmp; clear Ttmp;

h = m_O2*h_O2+m_N2*h_N2;
