function Psat = get_sat_partial_pressure(T)
global T_wi_max T_wi_min

% get expression for partial pressure of water vapor at saturation
Tc = 647.25; % critical temperature of water, K

if T<T_wi_min
    Psat = 611*exp(22.49-6142/T); % partial pressure of water at saturation
elseif (T<T_wi_max) 
    % in the temperature range (T_wi_min<T<T_wi_max) where liquid & ice coexist,
    % the partial pressure is weighted between the partial pressure of ice and of liquid.
    x = (T-T_wi_min)/(T_wi_max-T_wi_min);
    Psat = (1-x)*611*exp(22.49-6142/T)+x*1e5*exp(6.3573118-8858.843/T+607.56335/T^0.6);
elseif T<=314
    Psat = 1e5*exp(6.3573118-8858.843/T+607.56335/T^0.6);
elseif T>314 %&& T<Tc
    a = [-7.8889166 2.5514255 -6.716169 33.239495 -105.38479 174.35319 -148.39348 48.631602]; % fitting coefficients
    Psat = 2.2093e7*exp(Tc/T*sum(a.*(1-T/Tc).^(((1:8)+1)/2)));
    Psat = real(Psat);
end
  
end