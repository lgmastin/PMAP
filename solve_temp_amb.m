function Tamb = solve_temp_amb(Z)
global Z_trop H_trop dTdz_trop dTdz_strat Tamb_vent Z_vent

if Z < Z_trop % If we're still in the troposphere
    Tamb = Tamb_vent+dTdz_trop*(Z-Z_vent);
elseif Z < Z_trop+H_trop % in isothermal layer at tropopause
    Tamb = Tamb_vent+dTdz_trop*(Z_trop-Z_vent);
else % in stratosphere
    Tamb = Tamb_vent+dTdz_trop*(Z_trop-Z_vent)+dTdz_strat*((Z+Z_vent)-(Z_trop+H_trop));
end

end