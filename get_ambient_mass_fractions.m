function [m_v_amb,m_a_amb] = get_ambient_mass_fractions(Tamb,Pamb,RHamb)
global mw_air mw_H2O

P_H2O_sat = get_sat_partial_pressure(Tamb);

ws_amb = mw_H2O/mw_air*P_H2O_sat/(Pamb-P_H2O_sat); % mass ratio of water to dry air at saturation
w_amb = RHamb/100*ws_amb; % mass ratio of water vapor in humid air

m_v_amb = w_amb/(1+w_amb); % mass fraction of water in ambient
m_a_amb = 1-m_v_amb; % mass fraction of air in ambient
    
end