% Calculate Winding Conduction Loss for  4:1/2 SPVIRT Stackup Configurations
%
% Author:  Julia Estrin
% Date:    04-06-2026
%
% Description:
%   Computes the AC conduction losses in the primary and secondary windings
%   of a planar VIRT transformer using Dowell's method for various layer
%   stackup configurations. Supports 3-layer through 8-layer configurations,
%   both interleaved and non-interleaved.
%
%   Dowell coefficients:
%     F(delta) = [sinh(2d) + sin(2d)] / [cosh(2d) - cos(2d)]
%     G(delta) = [sinh(d)  - sin(d) ] / [cosh(d)  + cos(d) ]
%   where d = t_cu / skin_depth (normalized copper thickness).
%
%   Supported stackup configurations:
%     '3layer'             - 3-layer: P-P-S
%
% Inputs:
%   I         - Peak winding current [A]
%   f         - Switching frequency [Hz]
%   Rdc_pri   - Primary winding DC resistance [Ohm]
%   Rdc_sec   - Secondary winding DC resistance [Ohm]
%   t_cu_pri  - Primary copper trace thickness [m]
%   t_cu_sec  - Secondary copper trace thickness [m]
%   rho_cu    - Electrical resistivity of copper [Ohm·m]
%   u0        - Permeability of free space [H/m]
%   stackup   - Stackup configuration string (see above)
%
% Outputs:
%   P_pri     - Primary winding conduction loss [W]
%   P_sec     - Secondary winding conduction loss [W]
%
% Usage Example:
%   [P_pri, P_sec] = calculate_Pcond(10, 1e6, 0.05, 0.02, 70e-6, 70e-6, ...
%                                    1.72e-8, 4*pi*1e-7, '3layer');

function [P_pri, P_sec] = calculate_Pcond_SPVIRT(I, f, Rdc_pri, Rdc_sec, t_cu_pri, t_cu_sec, rho_cu, u0, stackup)

    %% Skin Depth and Normalized Copper Thickness
    skin_depth  = sqrt(rho_cu / (pi * f * u0));
    delta_pri   = t_cu_pri / skin_depth;
    delta_sec   = t_cu_sec / skin_depth;

    %% Dowell Coefficients
    F_pri = (sinh(2*delta_pri) + sin(2*delta_pri)) / (cosh(2*delta_pri) - cos(2*delta_pri));
    G_pri = (sinh(delta_pri)   - sin(delta_pri))   / (cosh(delta_pri)   + cos(delta_pri));
    
    F_sec = (sinh(2*delta_sec) + sin(2*delta_sec)) / (cosh(2*delta_sec) - cos(2*delta_sec));
    G_sec = (sinh(delta_sec)   - sin(delta_sec))   / (cosh(delta_sec)   + cos(delta_sec));

    %% Current Squared (Precompute)
    I2 = I^2;

    %% Calculate Losses Based on Stackup Configuration
    switch lower(stackup)
        
        case '3layer'
            % Configuration: P-P-S
            % 2 x factor for two parallel phases
           
            P_pri = 2 * 0.5 * delta_pri * Rdc_pri * (2 * I2 * F_pri + 4 * I2 * G_pri);
            P_sec = 2 * 0.5 * delta_sec * Rdc_sec * (4 * I2 *F_sec);
        
        
        
        otherwise
            error('Invalid stackup configuration: %s\nValid options: 3layer, 5layer, 6layer, 6layer_interleaved, 7layer_interleaved, 8layer_interleaved', stackup);
    end

end