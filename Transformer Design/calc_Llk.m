function Llk = calc_Llk(design_params, results)
t_pri    = design_params.t_cu_pri;
t_sec    = design_params.t_cu_sec;
t_shield = t_pri;

% Thickness of pre-preg 2116 (RC58% at 70% Cu fill) [m]
h_prepreg = 0.109e-3;
h1  = h_prepreg;
h3  = h_prepreg;
h5  = h_prepreg;
h7  = h_prepreg;
h9  = h_prepreg;
h11 = h_prepreg;

% Thickness of core [m]
h_core = 0.130e-3;
h2  = h_core;
h4  = h_core;
h6  = h_core;
h8  = h_core;
h10 = h_core;

w_wind = results.w_winding;
r_wind = (2*results.b + results.w_winding) / 2;  % average winding radius
l_wind = 2 * pi * r_wind;
u0     = design_params.u0;

Llk = 2 * u0 * l_wind / w_wind * ( ...
    16/3  * t_sec                          ...
    + 16  * (h10 + h11 + t_shield)        ...
    + 1/3 * t_pri + 1  * h7               ...
    + 7/3 * t_pri + 4  * h8               ...
    + 19/3* t_pri + 9  * h9               ...
    + 37/3* t_pri);
end