%% LLC Resonant Inductor Initial Design 
% Using core-area product and turn optimization from KPVS 
% JEstrin 05-11-2026

clear; clc;


iL_pk   = 13.79;        % [A]   Peak inductor current
Lr      = 4e-6;         % [H]   Resonant inductance
f       = 1e6;          % [Hz]  Switching frequency

B0      = 100e-3;       % [T]   Peak flux density (Ferrite 80 @ 1 MHz)
J0      = 2e6;          % [A/m²] Current density
ku      = 0.35;         % [-]   Winding fill factor

vL_pk   = 2*pi * f * Lr * iL_pk;   % [V]  Peak inductor voltage (V = L·dI/dt)
S       = 0.5 * vL_pk * iL_pk;     % [VA] Apparent power

% Calculate core-area product 
% Ac·Wa = S / (π·f·B0·J0·ku)   [m⁴]  →  convert to mm⁴
AcWa_m4  = S / (pi * f * B0 * J0 * ku);     % [m⁴]
AcWa_mm4 = AcWa_m4 * 1e12;                  % [mm⁴]

% Select 80 PLANAR E CORE22/6/16-80

Fr = 5; 

Pw = 