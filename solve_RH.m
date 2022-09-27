function RH = solve_RH(T,P,m_w,m_a)
global mw_H2O mw_air

if m_a <= 0
    RH = 0;
    return
end

P_H2O_sat = get_sat_partial_pressure(T);
P_H2O = m_w/mw_H2O/(m_w/mw_H2O+m_a/mw_air)*P; % partial pressure of water at plume conditions
RH = P_H2O/P_H2O_sat;

end