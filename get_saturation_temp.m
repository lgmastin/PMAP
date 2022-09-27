function Tsat = get_saturation_temp(P)
global R mw_H2O
% function that gives saturation temperature of water at the given partial pressure

Tsat = 182.82*P^0.0611; %first guess, using best-fit curve
Psat = get_sat_partial_pressure(Tsat);
while abs(Psat-P)>=10 % Adjust until we get within 10 Pascals
    Tsat = Tsat+(P-Psat)*R/mw_H2O*Tsat^2/((h_v(Tsat)-h_l(Tsat))*P);
    Psat = get_sat_partial_pressure(Tsat);
end

end