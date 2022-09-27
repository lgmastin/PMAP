function [T_next,m_v,m_l,m_i] = solve_temp(m_m,m_a,m_w,h,P)
global mw_H2O mw_air h_m T_wi_min T_wi_max hfl hfi hs Tsat
% subroutine that finds the mixture temperature given a mixture enthalpy and composition.

%SET DEFAULT VALUES
Cpa_avg = 1060;                 % average value of cp_air between 270 and 1000 K
Cpwv_avg = 2150;                % approximate average specific heat of water vapor between 100 c and 900 C
Cpwl_avg = 4190;                % Approximate average water specific heat
Cpwi_avg = 1850;                % average cp for ice
Cp_m = 1000;

h_tol = 1e-3; % allowable error in h when iterating temperature

if h>hs && h>hfl % if we're above water saturation, and above freezing
    m_v = m_w; m_l = 0; m_i = 0;
    Cp_mix = m_m*Cp_m+m_a*Cpa_avg+m_v*Cpwv_avg;
    % Estimate temperature, based on average specific heats.
    if hs>hfl
        T = Tsat+(h-hs)/Cp_mix;
    else
        T = T_wi_max+(h-hfl)/Cp_mix;
    end
    h_tmp = m_m*h_m(T)+m_v*h_v(T)+m_a*h_a(T);
    while abs(h_tmp-h)/h > h_tol
         %Iterate on final solution
         T = T+(h-h_tmp)/Cp_mix;
         h_tmp = m_m*h_m(T)+m_v*h_v(T)+m_a*h_a(T);
    end
elseif h>hfl && Tsat>T_wi_max % if we're within the saturated regime but above the freezing regime
    i = 1; m_i = 0;
    % Take a first stab at temperature
    T = T_wi_max+(Tsat-T_wi_max)*(h-hfl)/(hs-hfl);
    Psat = get_sat_partial_pressure(T);
    w_s = mw_H2O/mw_air*Psat/(P-Psat);
    if m_w >= w_s*m_a % we're above saturation
        m_v = w_s*m_a;
        m_l = m_w-m_v;
    else
        m_v = m_w;
        m_l = 0;
    end
    h_hi = hs;
    T_hi = Tsat;
    h_lo = hfl;
    T_lo = T_wi_max;
    if h_hi<h_lo 
        error('Problem with enthalpy calculations.')
    end
    h_tmp = m_m*h_m(T)+m_v*h_v(T)+m_l*h_l(T)+m_a*h_a(T);
    while abs(h_tmp-h)/h > h_tol % Iterate on final solution
        if h_tmp>h
            h_hi = h_tmp;
            T_hi = T;
        else
            h_lo = h_tmp;
            T_lo = T;
        end
        if h_hi<h_lo
            error('Problem with enthalpy ')
        end
        T = T_lo+(T_hi-T_lo)*(h-h_lo)/(h_hi-h_lo);
        if T<T_wi_max
            error('Problem with enthalpy ')
        end
        Psat = get_sat_partial_pressure(T);
        w_s = mw_H2O/mw_air*Psat/(P-Psat); % recalculate m_v
        if m_w >= w_s*m_a % we're above saturation
            m_v = w_s*m_a;
            m_l = m_w-m_v;
        else
            m_v = m_w;
            m_l = 0;
        end
        h_tmp = m_m*h_m(T)+m_v*h_v(T)+m_l*h_l(T)+m_a*h_a(T);
        i = i+1;
        %if i > 100 
        %    return
        %end
    end
elseif h>hfi % we're at temperature where ice and water may coexist
    i = 1;
    % Take a first stab at temperature
    T = T_wi_min+(T_wi_max-T_wi_min)*(h-hfi)/(hfl-hfi);
    Psat = get_sat_partial_pressure(T);
    w_s = mw_H2O/mw_air*Psat/(P-Psat);
    if w_s<0
        w_s = m_w/m_a; % w_s<0 only very high, perhaps at z>40 km
    end
    if m_w >= w_s*m_a % we're above saturation
        m_v = w_s*m_a;
        m_l = (m_w-m_v)*(T-T_wi_min)/(T_wi_max-T_wi_min);
        if m_l>1 
            disp('in FindT.  m_l>1.  m_l=',m_l)
            disp('pnow=',P,', psat(Tmix)=',psat(T))
            disp('m_a=',m_a,', w_s=',w_s)
            disp('m_w=',m_w,', m_v=',m_v)
            disp('Tmix=',T,', T_ice=',T_wi_min)
            disp('T_Coldwater=',T_wi_max)
            return
        end
        m_i = m_w-m_v-m_l;
    else
        m_v = m_w;
        m_l = 0;
        m_i = 0;
    end
    h_hi = hfl;
    T_hi = T_wi_max;
    h_lo = hfi;
    T_lo = T_wi_min;
    if h_hi<h_lo 
        error('Problem with enthalpy calculations.')
    end
    h_tmp = m_m*h_m(T)+m_v*h_v(T)+m_l*h_l(T)+m_i*h_i(T)+m_a*h_a(T);
    % Iterate on final solution
    T_prev = Inf;
    while abs(h_tmp-h)/h>h_tol && (T_hi - T_lo)>0.25 && abs(T_prev-T)>0.1
        if h_tmp>h 
            h_hi = h_tmp;
            T_hi = T;
        else
            h_lo = h_tmp;
            T_lo = T;
        end
        if h_hi<h_lo 
            error('Problem with enthalpy calculations.')
        end
        T_prev = T;
        T = T_lo+(T_hi-T_lo)*(h-h_lo)/(h_hi-h_lo);
        if T<T_wi_min
            error('Problem with enthalpy calculations.')
        end
        Psat = get_sat_partial_pressure(T);
        w_s = mw_H2O/mw_air*Psat/(P-Psat); % recalculate m_v
        if m_w >= w_s*m_a % we're above saturation
            m_v = w_s*m_a;
            m_l = (m_w-m_v)*(T-T_wi_min)/(T_wi_max-T_wi_min);
            m_i = m_w-m_v-m_l;
        else
            m_v = m_w;
            m_l = 0;
            m_i = 0;
        end
        h_tmp = m_m*h_m(T)+m_v*h_v(T)+m_l*h_l(T)+m_i*h_i(T)+m_a*h_a(T);
        i = i+1;
    end
else % if we're below pure ice temperature
    i = 1;
    m_l = 0;
    h_hi = hfi;
    h_lo = m_m*h_m(200)+m_w*h_i(200)+m_a*h_a(200);
    T_hi = T_wi_min;
    T_lo = 200;
    T = 200+(T_hi-T_lo)*(h-h_lo)/(h_hi-h_lo); % Take a first stab at temperature
    Psat = get_sat_partial_pressure(T);
    w_s = mw_H2O/mw_air*Psat/(P-Psat);
    if m_w >= w_s*m_a % we're above saturation
        m_v = w_s*m_a;
        m_i = m_w-m_v;
    else
        m_v = m_w;
        m_i = 0;
    end
    h_tmp = m_m*h_m(T)+m_v*h_v(T)+m_i*h_i(T)+m_a*h_a(T);
    while abs(h_tmp-h)/h>h_tol && (T_hi-T_lo)>0.25 % Iterate on final solution
        if h_tmp>h
            T_hi = T;
            h_hi = h_tmp;
        else
            T_lo = T;
            h_lo = h_tmp;
        end
        T = T_lo+(T_hi-T_lo)*(h-h_lo)/(h_hi-h_lo);
        Psat = get_sat_partial_pressure(T);
        w_s = mw_H2O/mw_air*Psat/(P - Psat); % recalculate m_v
        if m_w >= w_s*m_a % we're above saturation
            m_v = w_s*m_a;
            m_i = m_w-m_v;
        else
            m_v = m_w;
            m_i = 0;
        end
        h_tmp = m_m*h_m(T)+m_v*h_v(T)+m_i*h_i(T)+m_a*h_a(T);
%         if (i > 20)
%             return
%         end
    end
end 

T_next = T;

end