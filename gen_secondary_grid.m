function [M_grid,rho_grid,ETA] = gen_secondary_grid
global Mbin rho_bin Nm Nrho Nc Df

% create 2D staggered (triangular) grid based on specified M and rho spacing
M_grid = Mbin;
M_centers = (Mbin(1:end-1)+Mbin(2:end))/2;
for i = 2:Nrho
    M_grid = [M_grid M_centers];
    M_grid = [M_grid Mbin];
end
rho_grid = rho_bin(1)*ones(1,Nm);
rho_centers = (rho_bin(1:end-1)+rho_bin(2:end))/2;
for i = 2:Nrho
    rho_grid = [rho_grid rho_centers(i-1)*ones(1,Nm-1)];
    rho_grid = [rho_grid rho_bin(i)*ones(1,Nm)];
end

% create voronoi tesselation and connectivity list of grid
DT = delaunayTriangulation(M_grid',rho_grid');
C = DT.ConnectivityList;

%%
% create (Nc x Nc) lookup matrix to store the voronoi index for every particle class combo
disp('Generating look up table for secondary grid')
Mmax = max(M_grid);
rho_min = min(rho_grid);

tic
lookup = zeros(Nc);
for j = 1:Nc
    for k = 1:j
        M12 = M_grid(j)+M_grid(k);
        rho12 = M12*((M_grid(j)/rho_grid(j))^(Df/3)+(M_grid(k)/rho_grid(k))^(Df/3))^(-3/Df);
        if M12<Mmax && rho12>rho_min
            lookup(j,k) = pointLocation(DT,M12,rho12);
            lookup(k,j) = lookup(j,k);
        end
    end
end
elapsed = toc;
disp(['Generated look-up table in ' num2str(elapsed) ' s'])

%%
% check if weights matrix for the current settings is already stored
weights_matrix = 'input\weights.mat';
if isfile(weights_matrix)
    load(weights_matrix)
    if D_stored==Df && sum(Mbin_stored)==sum(Mbin) && sum(rho_bin_stored)==sum(rho_bin)
        return
    end
    clear D_stored Mbin_stored rho_bin_stored
end

%%
% generate (Nc x Nc x Nc) weighting matrix
% mass fraction of particle i/j combo that ends up in pivot k
disp('Generating weighting matrix for pivots')

tic
ETA = zeros(Nc,Nc,Nc);
for i = 1:Nc
    for j = 1:Nc
        for k = 1:j
            eta = 0;
            M12 = M_grid(j)+M_grid(k);
            rho12 = M12*((M_grid(j)/rho_grid(j))^(Df/3)+(M_grid(k)/rho_grid(k))^(Df/3))^(-3/Df);
            if M12<Mbin(end) && rho12>rho_bin(end)
                ID = lookup(j,k);
                if any(C(ID,:)==i)
                    eta = get_eta_internal(i,C(ID,:),M12,rho12,M_grid,rho_grid);
                end
            elseif M12>=Mbin(end) && rho12>rho_bin(end)
                rho_ind = find((rho12-rho_bin)>0);
                rho_ind = [rho_ind(1) rho_ind(1)-1];
                if M_grid(i)==Mbin(end) && any(rho_grid(i)==rho_bin(rho_ind))
                    eta = get_eta_Medge(i,rho_ind,rho12,rho_grid);
                end
            elseif rho12<=rho_bin(end) && M12<Mbin(end)
                M_ind = find((Mbin-M12)>0);
                M_ind = [M_ind(1) M_ind(1)-1];
                if rho_grid(i)==rho_bin(end) && any(M_grid(i)==Mbin(M_ind))
                    eta = get_eta_rho_edge(i,M_ind,M12,M_grid);
                end
            elseif M12>=Mbin(end) && rho12<=rho_bin(end)
                if M_grid(i)==Mbin(end) && rho_grid(i)==rho_bin(end)
                    eta = M12/M_grid(i);
                end
            end
            ETA(i,j,k) = eta;
            ETA(i,k,j) = eta;
        end
    end
end
timer = toc;

disp(['Generated weighting matrix in ' num2str(timer) ' s'])
D_stored = Df; Mbin_stored = Mbin; rho_bin_stored = rho_bin;
save(weights_matrix,'D_stored','Mbin_stored','rho_bin_stored','ETA')


function eta = get_eta_internal(i,C,M,rho,M_grid,rho_grid)

C(C==i) = [];

b = ones(3,1);
A = [rho_grid(i) rho_grid(C(1)) rho_grid(C(2));
     M_grid(i)   M_grid(C(1))   M_grid(C(2));
     1           1              1        ];
A(1,:) = A(1,:)/rho;
A(2,:) = A(2,:)/M;

w = A\b;
eta = w(1);

end

% Helper functions for weighting at the edges of the grid
function eta = get_eta_Medge(i,rho_ind,rho,rho_grid)
% global rho_bin

i = find(rho_grid(i)==rho_bin(rho_ind));
eta = (rho-rho_bin(rho_ind(3-i)))/(rho_bin(rho_ind(i))-rho_bin(rho_ind(3-i)));

end

function eta = get_eta_rho_edge(i,M_ind,M,M_grid)
% global Mbin

i = find(M_grid(i)==Mbin(M_ind));
eta = (M-Mbin(M_ind(3-i)))/(Mbin(M_ind(i))-Mbin(M_ind(3-i)));

end

end