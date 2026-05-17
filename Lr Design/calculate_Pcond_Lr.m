% Calculate Winding Conduction for Multitrack Resonant Inductor 
%
% Author:  Julia Estrin
% Date:    05-11-2026
%
% Description:
%   Computes the AC conduction losses in planar inducuctor using Dowell's 
%   method 

%   Dowell coefficients:
%     F(delta) = [sinh(2d) + sin(2d)] / [cosh(2d) - cos(2d)]
%     G(delta) = [sinh(d)  - sin(d) ] / [cosh(d)  + cos(d) ]
%   where d = t_cu / skin_depth (normalized copper thickness).
%

function [P_cond] = calculate_Pcond_Lr(I, f, Rdc, N, t_cu, rho_cu, u0)

%% Skin Depth and Normalized Copper Thickness
skin_depth  = sqrt(rho_cu / (pi * f * u0));
delta   = t_cu / skin_depth;

%% Dowell Coefficients
F_delta = (sinh(2*delta) + sin(2*delta)) / (cosh(2*delta) - cos(2*delta));
G_delta = (sinh(delta)   - sin(delta))   / (cosh(delta)   + cos(delta));

    P_cond = 0; 
    for n = 1:N 
        
        P_cond = 1/2* Rdc * delta * ((n*I - (n-1)*I)^2 * F_delta ...
            + 2 * n*I * (n-1)*I * G_delta); 
    
    end 

end