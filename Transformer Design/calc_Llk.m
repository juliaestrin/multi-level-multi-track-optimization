function Llk = calc_Llk(design_params, results)
    t_ins = 0.346e-3; 
    t_pri = design_params.t_cu_pri; 
    t_sec = design_params.t_cu_sec; 

    w_wind = results.w_winding; 
    r_wind = (2*results.b + results.w_winding)/2;
    l_wind = 2*pi*r_wind; 

    u0 = design_params.u0; 

    Llk = u0 * l_wind * w_wind * (8 *t_pri/(3 * w_wind^2) + 56 * t_pri / (3 * w_wind^2) ...
        + 8 * t_ins/w_wind^2 + 32 * t_ins/w_wind^2 ...
        + 16 * t_sec/(3*w_wind^2));

end 