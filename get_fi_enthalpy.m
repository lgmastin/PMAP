function hfi = get_fi_enthalpy(m_w,m_a,m_m,P)
% find enthalpy at saturation and all excess water is ice
global mw_H2O mw_air h_m T_wi_min

P_H2O_sat = get_sat_partial_pressure(T_wi_min); % partial pressure of water at fully frozen temp

if P>P_H2O_sat % if the ambient pressure exceeds the boiling pressure at T_ice
    ws_fi = mw_H2O/mw_air*P_H2O_sat/(P-P_H2O_sat);
    m_v_fi = ws_fi/(1+ws_fi);
    if m_w>m_v_fi % if there's enough water vapor to saturate the plume
        hfi = m_m*h_m(T_wi_min)+m_a*h_a(T_wi_min)+m_v_fi*h_v(T_wi_min)+(m_w-m_v_fi)*h_i(T_wi_min);
    else % if air is non water-saturated at T_ice
        m_v_fi = m_w;
        hfi = m_m*h_m(T_wi_min)+m_a*h_a(T_wi_min)+m_v_fi*h_v(T_wi_min);
    end
else                                
    % if pnow is less than the boiling pressure at T_ice
    m_v_fi = m_w;
    hfi = m_m*h_m(T_wi_min)+m_a*h_a(T_wi_min)+m_v_fi*h_v(T_wi_min);
end
          
end