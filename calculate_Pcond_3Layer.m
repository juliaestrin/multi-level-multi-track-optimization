% Calculate Winding Conduction Loss for 3-Layer Planar VIRT Windings
%
% Author:  Julia Estrin
% Date:    02-12-2026
%
% Description:
%   Computes the AC conduction losses in the primary and secondary windings
%   of a 3-layer planar VIRT transformer using Dowell's method. The skin
%   and proximity effect losses are captured via the Dowell coefficients
%   F(delta) and G(delta), which are functions of the normalized copper
%   thickness (ratio of physical thickness to skin depth).
%
%   Dowell coefficients:
%     F(delta) = [sinh(2d) + sin(2d)] / [cosh(2d) - cos(2d)]
%     G(delta) = [sinh(d)  - sin(d) ] / [cosh(d)  + cos(d) ]
%   where d = t_cu / skin_depth (normalized copper thickness).
%
%   Loss expressions:
%     P_pri = 4 * delta_pri * Rdc_pri * I^2 * F(delta_pri)
%     P_sec = (1/2) * Rdc_sec * delta_sec * 8 * I^2 * (2*F(delta_sec) - G(delta_sec))
%
% Inputs:
%   I         - RMS winding current [A]
%   f         - Switching frequency [Hz]
%   Rdc_pri   - Primary winding DC resistance [Ohm]
%   Rdc_sec   - Secondary winding DC resistance [Ohm]
%   t_cu_pri  - Primary copper trace thickness [m]
%   t_cu_sec  - Secondary copper trace thickness [m]
%   rho_cu    - Electrical resistivity of copper [Ohm·m]  (typically ~1.72e-8)
%   u0        - Permeability of free space [H/m]          (4*pi*1e-7)
%
% Outputs:
%   P_pri     - Primary winding conduction loss [W]
%   P_sec     - Secondary winding conduction loss [W]
%

function [P_pri, P_sec] = calculate_Pcond_3Layer(I, f, Rdc_pri, Rdc_sec, t_cu_pri, t_cu_sec, rho_cu, u0)

    % --- Skin depth and normalized copper thickness ---
    skin_depth  = sqrt(rho_cu / (pi * f * u0));
    delta_pri   = t_cu_pri / skin_depth;
    delta_sec   = t_cu_sec / skin_depth;

    % --- Dowell coefficients F(delta) and G(delta) ---
    F_pri = (sinh(2*delta_pri) + sin(2*delta_pri)) / (cosh(2*delta_pri) - cos(2*delta_pri));
    G_pri = (sinh(delta_pri)   - sin(delta_pri))   / (cosh(delta_pri)   + cos(delta_pri));

    F_sec = (sinh(2*delta_sec) + sin(2*delta_sec)) / (cosh(2*delta_sec) - cos(2*delta_sec));
    G_sec = (sinh(delta_sec)   - sin(delta_sec))   / (cosh(delta_sec)   + cos(delta_sec));

    % --- Conduction losses ---
    P_pri   = 4   * delta_pri * Rdc_pri * I^2 * F_pri;
    P_sec   = 8   * delta_sec * Rdc_sec * I^2 * (2*F_sec - G_sec);

end