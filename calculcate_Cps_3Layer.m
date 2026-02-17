function Cps = calculcate_Cps_3Layer(l_winding, w_winding, w_core)
    e0 = 8.854e-12; % peritivity of free space
    er = 4.2; % permitivity of FR4
    d_ij = 65e-6; % [m] distance between pcb layers 
    
    A_ij = l_winding*w_winding - w_core^2; % [m^2] winding area 

    C_ij = er*e0*A_ij/d_ij; 

    Cps = C_ij * 10/9;
end
