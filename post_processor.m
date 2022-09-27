close all

% specify file to read here if running stand-alone
if exist('savefilename','var')~=1
    load(savefilename)
end

plt_clr = 'r'; % plot color specification

Nc = numel(output.rho_p);
Nm = sum(output.rho_p==max(output.rho_p));
dp = (6/pi*output.Mp./output.rho_p).^(1/3);
phi = -log2(dp/1e-3);
phi_max = max(phi);
phi_min = min(phi);

%% Plots
i = 1;

% check that all mass has been accounted for
figure(i)
Maccounting = (sum(output.Mp.*pi.*output.r.^2.*output.u.*output.Nd,2)+cumsum([0; output.Mfall]));
fprintf('Premature fall-out of tephra is %.1f%% by mass\n',sum(output.Mfall)/Maccounting(1)*100)
plot(Maccounting/Maccounting(1),output.Z/1e3,plt_clr)
xlabel('Mass balance (%)')
ylabel('$z$ (km)','interpreter','latex')
scale = axis; axis([scale(1) scale(2) 0 scale(4)])
i = i+1;

% plot ascent velocity
figure(i)
plot(output.u,output.Z/1e3,plt_clr); hold on
xlabel('$u$ (m/s)','interpreter','latex')
ylabel('$z$ (km)','interpreter','latex')
scale = axis; axis([scale(1) scale(2) 0 scale(4)])
i = i+1;

% plot plume radius
figure(i)
plot(output.r/1e3,output.Z/1e3,plt_clr); hold on
xlabel('$r$ (km)','interpreter','latex')
ylabel('$z$ (km)','interpreter','latex')
scale = axis; axis([scale(1) scale(2) 0 scale(4)])
i = i+1;

% plot mass fractions
figure(i)
plot(output.m_w,output.Z/1e3,'k-'); hold on
plot(output.m_a,output.Z/1e3,'g-')
plot(output.m_m,output.Z/1e3,'r-')
xlabel('mass fraction')
ylabel('$z$ (km)','interpreter','latex')
legend('water','air','magma'); legend boxoff
scale = axis; axis([scale(1) scale(2) 0 scale(4)])
i = i+1;

% plot pressure
figure(i)
plot(output.P,output.Z/1e3,plt_clr); hold on
xlabel('$P$ (Pa)','interpreter','latex')
ylabel('$z$ (km)','interpreter','latex')
scale = axis; axis([scale(1) scale(2) 0 scale(4)])
i = i+1;

% plot temperature
figure(i)
plot(output.T,output.Z/1e3,plt_clr); hold on
plot(output.Tamb,output.Z/1e3,[plt_clr '--'])
xlabel('$T$ (K)','interpreter','latex')
ylabel('$z$ (km)','interpreter','latex')
legend('Plume','Ambient','Location','NorthEast'); legend boxoff
scale = axis; axis([scale(1) scale(2) 0 scale(4)])
i = i+1;

% plot density
figure(i)
plot(output.rho,output.Z/1e3,plt_clr); hold on
plot(output.rho_amb,output.Z/1e3,[plt_clr '--'])
xlabel('$\rho$ (kg/m$^3$)','interpreter','latex')
ylabel('$z$ (km)','interpreter','latex')
legend('Plume','Ambient','Location','NorthEast'); legend boxoff
scale = axis; axis([scale(1) scale(2) 0 scale(4)])
i = i+1;

% plot hydrous mass fraction
figure(i)
semilogx(output.m_v,output.Z/1e3,'k-'); hold on
plot(output.m_l,output.Z/1e3,'g-'); 
plot(output.m_i,output.Z/1e3,'b-'); 
xlabel('hydrous mass fraction')
ylabel('$z$ (km)','interpreter','latex')
legend('vapor','liquid','ice'); legend boxoff
scale = axis; axis([scale(1) scale(2) 0 scale(4)])
i = i+1;

 % plot grain size distribution
figure(i)
Gm_f_top = zeros(1,Nm);
Gm_f_tot = zeros(1,Nm);
for j = 1:Nc
    [~,ind] = min(abs(dp(j)-dp(1:Nm)));
    Gm_f_top(ind) = Gm_f_top(ind)+output.Mp(j)*output.n_detrain(end,j);
    Gm_f_tot(ind) = Gm_f_tot(ind)+output.Mp(j)*sum(output.n_detrain(:,j));
end
Gm_i = output.Nd(1,1:Nm).*dp(1:Nm).^3/sum(output.Nd(1,1:Nm).*dp(1:Nm).^3);
Gm_f_top = Gm_f_top/sum(Gm_f_top);
Gm_f_tot = Gm_f_tot/sum(Gm_f_tot);
plot(phi(1:Nm),Gm_i,'k-o'); hold on
plot(phi(1:Nm),Gm_f_tot,'k--d');
xlabel('$\phi$','interpreter','latex')
ylabel('mass fraction','interpreter','latex')
legend('Initial GSD','AGSD','Location','Northwest'); legend boxoff
scale = axis; axis([phi_min phi_max 0 scale(4)])
i = i+1;

if input.Df~=3
    figure(i); hold on
    tri = delaunay(output.Mp,output.rho_p);
    Gm = output.Mp.*sum(output.n_detrain,1);
    Gm = Gm/sum(Gm);
    for j = 1:Nc
        plot3(output.Mp(j),output.rho_p(j),Gm(j),'ko','MarkerFaceColor','k')
        plot3(output.Mp(j)*ones(2,1),output.rho_p(j)*ones(2,1),Gm(j)*[0 1],'k-','LineWidth',2)
    end
    ts = trisurf(tri,output.Mp,output.rho_p,Gm); colormap hot
    ts.FaceAlpha = 0.75;
    set(gca, 'XScale', 'log')
    view(22.078,63.566)
    xlabel('$M_\mathrm{p}$ (kg)','interpreter','latex')
    yl = ylabel('$\rho_\mathrm{p}$ (kg/m$^3$)','interpreter','latex');
    zlabel('mass fraction')
    shading interp

    DT = delaunayTriangulation(output.Mp',output.rho_p');
    triplot(DT,'-','Color',0*[1 1 1],'LineWidth',1)

    % plot parent GSD
    load('input\Redoubt_Event5_GSD.mat'); Gm = GSD; clear GSD
    Gm = Gm/sum(Gm);
    plot3(output.Mp(1:Nm),output.rho_p(1:Nm),Gm,'-o','Color',0.5*[1 1 1],'MarkerFaceColor',0.5*[1 1 1],'LineWidth',2)
end

% plotfixer