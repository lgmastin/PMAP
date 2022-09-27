function hfl = get_fl_enthalpy(m_w,m_a,m_m,P)
% find enthalpy at saturation and all excess water is liquid
global mw_H2O mw_air h_m T_wi_max

P_H2O_sat = get_sat_partial_pressure(T_wi_max); % partial pressure of water at fully frozen temp

if P>P_H2O_sat % if the ambient pressure exceeds the boiling pressure at T_wi_max
    ws_fl = mw_H2O/mw_air*P_H2O_sat/(P-P_H2O_sat);
    m_v_fl = ws_fl/(1+ws_fl);
    if m_w>m_v_fl % if there s enough water vapor to saturate the plume
        hfl = m_m*h_m(T_wi_max)+m_a*h_a(T_wi_max)+m_v_fl*h_v(T_wi_max)+(m_w-m_v_fl)*h_l(T_wi_max);
    else %if air is non water-saturated at T_ice
        m_v_fl = m_w;
        hfl = m_m*h_m(T_wi_max)+m_a*h_a(T_wi_max)+m_v_fl*h_v(T_wi_max);
    end
else                                
    % if pnow is less than the boiling pressure at T_ice
    m_v_fl = m_w;
    hfl = m_m*h_m(T_wi_max)+m_a*h_a(T_wi_max)+m_v_fl*h_v(T_wi_max);
end
          
end