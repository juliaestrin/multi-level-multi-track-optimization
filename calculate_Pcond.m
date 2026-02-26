% Calculate Winding Conduction Loss for Various VIRT Stackup Configurations
%
% Author:  Julia Estrin
% Date:    02-12-2026
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
%     '3layer'             - 3-layer: P-S-P
%     '5layer'             - 5-layer: P-P-P-P-S
%     '5layer_interleaved' - 5-layer interleaved: P-P-S-P-P
%     '6layer'             - 6-layer non-interleaved: P-P-S-S-P-P
%     '6layer_interleaved' - 6-layer interleaved: P-S-P-P-S-P
%     '7layer_interleaved' - 7-layer interleaved: P-S-P-S-P-S-P
%     '8layer_interleaved' - 8-layer interleaved: P-S-P-S-P-S-P-S
%
% Inputs:
%   I         - RMS winding current [A]
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

function [P_pri, P_sec] = calculate_Pcond(I, f, Rdc_pri, Rdc_sec, t_cu_pri, t_cu_sec, rho_cu, u0, stackup)

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
            % Configuration: P-S-P
            % Primary outer layers (2 layers, each sees only skin effect)
            % Secondary center layer (sees proximity effect from both sides)
            P_pri = 2 * 0.5 * delta_pri * Rdc_pri * (4 * I2 * F_pri);
            P_sec = 0.5 * delta_sec * Rdc_sec * I2 * (16*F_sec - 8*G_sec);
        
        case '5layer'
            % Configuration: P-P-P-P-S
            P_pri = 0.5 * delta_pri * Rdc_pri * (4 * I2 * F_pri + 2 * (0 + 2 + 6 + 12) * G_pri);
            
            % Secondary: center layer with strong proximity from 4 primary layers
            P_sec = 0.5 * delta_sec * Rdc_sec * I2 * (16*F_sec);

        case '5layer_interleaved'
            % Configuration: P-P-S-P-P
            % Primary: 2 outer layers + 2 inner layers
            % Outer layers: skin effect only
            % Inner layers: skin + proximity effect (m=2 position)
            P_pri_outer = 2 * 0.5 * delta_pri * Rdc_pri * I2 * F_pri;
            P_pri_inner = 2 * 0.5 * delta_pri * Rdc_pri * I2 * (F_pri + 4*G_pri);
            P_pri = P_pri_outer + P_pri_inner;
            
            % Secondary: center layer with strong proximity from 4 primary layers
            P_sec = 0.5 * delta_sec * Rdc_sec * I2 * (16*F_sec - 8*G_sec);
        
        case '6layer'
            % Configuration: P-P-S-S-P-P (non-interleaved)
            % Primary: 2 outer + 2 inner layers
            P_pri_outer = 2 * 0.5 * delta_pri * Rdc_pri * I2 * F_pri;
            P_pri_inner = 2 * 0.5 * delta_pri * Rdc_pri * I2 * (F_pri + 4*G_pri);
            P_pri = P_pri_outer + P_pri_inner;
            
            % Secondary: 2 center layers, each sees proximity but reduced G term
            P_sec = 2 * 0.5 * delta_sec * Rdc_sec * I2 * (4*F_sec);
        
        case '6layer_interleaved'
            % Configuration: P-S-P-P-S-P (interleaved)
            % Primary: 2 outer + 2 inner, all see same field distribution
            P_pri_outer = 2 * 0.5 * delta_pri * Rdc_pri * I2 * F_pri;
            P_pri_inner = 2 * 0.5 * delta_pri * Rdc_pri * I2 * F_pri;
            P_pri = P_pri_outer + P_pri_inner;
            
            % Secondary: 2 layers interleaved, reduced proximity effect
            P_sec = 2 * 0.5 * delta_sec * Rdc_sec * I2 * (4*F_sec - 2*G_sec);
        
        case '7layer_interleaved'
            % Configuration: P-S-P-S-P-S-P (interleaved)
            % Primary: 2 outer + 2 inner layers
            P_pri_outer = 2 * 0.5 * delta_pri * Rdc_pri * I2 * F_pri;
            P_pri_inner = 2 * 0.5 * delta_pri * Rdc_pri * I2 * F_pri;
            P_pri = P_pri_outer + P_pri_inner;
            
            % Secondary: 3 layers interleaved
            % Note: Formula from document has 6*F - 2*G in the bracket
            P_sec = 0.5 * delta_sec * Rdc_sec * I2 * (6*F_sec - 2*G_sec);
        
        case '8layer_interleaved'
            % Configuration: P-S-P-S-P-S-P-S (fully interleaved)
            % Primary: 2 outer + 2 inner layers
            P_pri_outer = 2 * 0.5 * delta_pri * Rdc_pri * I2 * F_pri;
            P_pri_inner = 2 * 0.5 * delta_pri * Rdc_pri * I2 * F_pri;
            P_pri = P_pri_outer + P_pri_inner;
            
            % Secondary: 4 layers interleaved, minimal proximity effect
            P_sec = 2 * 0.5 * delta_sec * Rdc_sec * I2 * (2*F_sec);
        
        otherwise
            error('Invalid stackup configuration: %s\nValid options: 3layer, 5layer, 6layer, 6layer_interleaved, 7layer_interleaved, 8layer_interleaved', stackup);
    end

end