function [dydx] = derivs(x,yarr)
global R mw_H2O mw_air h_m g RHamb rho_m rho_l rho_i eps Tsat hs hfi hfl ...
    emmisivity sigma Chi up_t Nc metdata Zmet Tmet Mp rho_p rp dp Nm Df ...
    AggModel Ko ETA Mbin aggregation kB StickModel

Z = x;

% store variables from yarr
M   = yarr(1);
p   = yarr(2);
E   = yarr(3);
P   = yarr(4);
m_m = yarr(5);
m_a = yarr(6);
m_w = yarr(7);
n  = yarr(8:end);

if metdata
    Tamb = interp_custom(Zmet,Tmet,Z);
else
    Tamb = solve_temp_amb(Z); % ambient temperature
end
[m_v_amb,m_a_amb] = get_ambient_mass_fractions(Tamb,P,RHamb); % ambient composition
h_amb = m_a_amb*h_a(Tamb)+m_v_amb*h_v(Tamb); % ambient enthalpy
rho_amb = P/((m_v_amb/mw_H2O+m_a_amb/mw_air)*R*Tamb); % ambient density

u = p/M; % axial velocity
h = (E/M)-u^2/2-g*Z; % mixture enthalpy
Tsat = get_saturation_temp(m_w/mw_H2O/(m_w/mw_H2O+m_a/mw_air)*P); % saturation temperature
hs = m_w*h_v(Tsat)+m_a*h_a(Tsat)+m_m*h_m(Tsat); % saturation enthalpy
hfi = get_fi_enthalpy(m_w,m_a,m_m,P); % freezing enthalpy
hfl = get_fl_enthalpy(m_w,m_a,m_m,P); % condensed enthalpy
 
[T,m_v,m_l,m_i] = solve_temp(m_m,m_a,m_w,h,P); % mixture temperature
rho = 1/(m_l/rho_l+m_m/rho_m+m_i/rho_i+(m_a/mw_air+m_v/mw_H2O)*R*T/P); % mixture density
r = sqrt(M/(pi*u*rho)); % plume radius

Nd = n/(pi*r^2*u);
Mc = pi/6*dp.^3.*rho_p.*n;
dMcdz = -Chi*up_t/(r*u).*Mc;
dndz = -Chi*up_t/(r*u).*n;

% get water/ice film thickness
tliq = 0;
tice = 0;
if m_l+m_i>0
    % get starting value assuming thin film
    tliq = M*m_l/(4*pi*rho_l*sum(n.*rp.^2));
    tice = M*m_i/(4*pi*rho_i*sum(n.*rp.^2));

    % get accurate value using iterative solver
    func = @(t) 4/3*pi*rho_l*sum(n.*((rp+t).^3-rp.^3))-M*m_l;
    tliq = fzero(func,tliq);
    func = @(t) 4/3*pi*rho_i*sum(n.*((rp+tliq+t).^3-(rp+tliq).^3))-M*m_i;
    tice = fzero(func,tice);
end
Mfilm_liq = 4/3*pi*((rp+tliq).^3-rp.^3)*rho_l;
Mfilm_ice = 4/3*pi*((rp+tliq+tice).^3-(rp+tliq).^3)*rho_i;

% calculate derivatives
if rho>rho_amb
    u_e = sqrt(rho/rho_amb)*eps*u;
else
    u_e = eps*u;
end
dMdz = 2*pi*r*rho_amb*u_e+sum(dndz.*(Mp+Mfilm_liq+Mfilm_ice));
dydx(1) = dMdz; % dM/dz
dpdz = pi*r^2*g*(rho_amb-rho);  % dMom/dz
dydx(2) = dpdz;
dydx(3) = 2*pi*r*rho_amb*u_e*(h_amb+u_e^2/2+g*Z)+(h_m(T)+g*Z)*sum(dndz.*Mp)+ ...
    (h_l(T)+g*Z)*sum(dndz.*Mfilm_liq)+(h_i(T)+g*Z)*sum(dndz.*Mfilm_ice)- ...
    2*pi*r*emmisivity*sigma*(T^4-Tamb^4); % dE/dz
dydx(4) = -rho_amb*g; % dP/dz
dydx(5) = -sum(dndz.*Mp)/M-dMdz*m_m/M; % dm_m/dz
dydx(6) = 2*pi*r*rho_amb*u_e*m_a_amb/M-dMdz*m_a/M; % dm_a/dz
dydx(7) = (2*pi*r*rho_amb*u_e*m_v_amb-sum(dndz.*Mfilm_liq)-sum(dndz.*Mfilm_ice))/M-dMdz*m_w/M; % dm_w/dz

%%
smol = zeros(1,Nc);
dydx(8:8+Nc-1) = 0;
RH = solve_RH(T,P,m_w,m_a);
Magg_tot = 0;
Magg_DS = 0;
Magg_TI = 0;
if aggregation && m_l+m_i~=0
    % get properties of plume gas phase and calculate particle settling velocities
    mu = 1.458e-6*T.^1.5./(T+110.4); % dynamic viscosity, kg/m.s
    m_a_g = m_a/(m_a+m_v);
    m_v_g = m_v/(m_a+m_v);
    rho_g = P/(m_a_g/mw_air+m_v_g/mw_H2O)/R/T;
    nu = mu./rho_g; % kinematic viscosity, m^2/s
    epsilon = 0.017*u^3/r; % dissipation rate, m^2/s^3 (Panchapakesan and Lumley, 1993)
    eta = (nu^3/epsilon)^(1/4);
    tau_eta = sqrt(nu/epsilon);
    tau_p = up_t/g;
    St_eta = tau_p/tau_eta;
    
    % define bidisperse RDF
    phi = @(i,j) max([tau_p(i) tau_p(j)])/min([tau_p(i) tau_p(j)]);
    corr = @(i,j) 2.6*exp(-phi(i,j))+0.205*exp(-0.0206*phi(i,j))*(1+tanh(phi(i,j)-3))/2;
    c0 = (7.92*St_eta.^1.80)./(0.58+St_eta.^3.29);
    c1 = (0.61*St_eta.^0.88)./(0.33+St_eta.^2.38);
    gmono = @(i) c0(i)*(dp(i)/eta)^(-c1(i))*exp(-0.25*dp(i)/eta)+1;
    g_bi = @(i,j) 1+corr(i,j)*sqrt((gmono(i)-1)*(gmono(j)-1));
    
    % define collision kernel
    K = @(i,j) AggModel(1)*Ko+ ...
               AggModel(2)*(rp(i)+rp(j))^3*sqrt(8*pi/15*epsilon/nu)+ ...
               AggModel(3)*(rp(i)+rp(j))^3*g_bi(i,j)*sqrt(8*pi/15*epsilon/nu)+ ...
               AggModel(4)*pi*(rp(i)+rp(j))^2*abs(up_t(i)-up_t(j))+ ...
               AggModel(5)*2*kB*T/3/mu*(rp(i)+rp(j))^2/rp(i)/rp(j);
                    
    % define viscous Stokes number
    Ur = @(i,j) (AggModel(2)+AggModel(3))*(rp(i)+rp(j))*sqrt(2/15/pi*epsilon/nu)+ ...
                AggModel(4)*abs(up_t(i)-up_t(j))+ ...
                AggModel(5)*kB*T/3/pi/mu/rp(i)/rp(j);
    
    % define wet sticking efficiency
    if strcmp(StickModel,'Costa')
        St = @(i,j) 16/9*(rho_p(i)+rho_p(j))/2/(1e-6*1000)*Ur(i,j)*rp(i)*rp(j)/(rp(i)+rp(j)); % Stokes number
        alpha_liq = @(i,j) 1/(1+(St(i,j)/1.3)^0.8);
    elseif strcmp(StickModel,'Hoffman')
        St = @(i,j) 2/(3*pi)*Mp(i)*Mp(j)/(Mp(i)+Mp(j))/8.9e-4*Ur(i,j)/(dp(i)*dp(j)/(dp(i)+dp(j)))^2; % Stokes number
        alpha_liq = @(i,j) 1/(1+(St(i,j)/max(0,2*log(tliq/6.5e-9)))^0.9);
    end
    
    % define combined sticking efficiency
    alpha_ice = @(i,j) 0.09; % assume constant sticking efficiency for ice
    if m_l>0 || m_i>0 % saturated conditions with liquid water and/or ice 
        alpha = @(i,j) (m_l*alpha_liq(i,j)+m_i*alpha_ice(i,j))/(m_l+m_i);
    else
        alpha = @(i,j) 0;
    end
    
    % Smoluchowski equation
    if Df == 3  % coalescing aggregation
        for i = 1:Nm 
            bsum = 0;
            for j = 1:Nm
                for k = 1:Nm
                    Mkj = Mbin(j)+Mbin(k);
                    eta = get_eta_coal(i,Mkj);
                    if eta~=0
                        bsum = bsum+eta*alpha(j,k)*K(j,k)*n(j)*Nd(k);
                    end
                end
            end
            dsum = 0;
            for k = 1:Nm
                dsum = dsum+alpha(i,k)*K(i,k)*Nd(k);
            end
            smol(1,i) = 1/u*(1/2*bsum-n(i)*dsum);
        end
    else % fractal aggregation
        for i = 1:Nc 
            bsum = 0;
            for j = 1:Nc
                for k = 1:Nc
                    eta = ETA(i,j,k);
                    if eta~=0
                        bsum = bsum+eta*alpha(j,k)*K(j,k)*n(j)*Nd(k);
                    end
                end
            end
            dsum = 0;
            for k = 1:Nc
                dsum = dsum+alpha(i,k)*K(i,k)*Nd(k);
            end
            smol(1,i) = 1/u*(1/2*bsum-n(i)*dsum);
        end
    end
end
dydx(8:8+Nc-1) = smol-Chi*up_t/r/u.*n;

%%
function eta = get_eta_coal(i,Mkj)

if i<Nm && Mkj>=Mbin(i) && Mkj<=Mbin(i+1)
    eta = (Mbin(i+1)-Mkj)/(Mbin(i+1)-Mbin(i));
elseif i>1 && Mkj>=Mbin(i-1) && Mkj<=Mbin(i)
    eta = (Mkj-Mbin(i-1))/(Mbin(i)-Mbin(i-1));
elseif i==Nm && Mkj>=Mbin(i)
    eta = 1;
else
    eta = 0;
end

end

end
