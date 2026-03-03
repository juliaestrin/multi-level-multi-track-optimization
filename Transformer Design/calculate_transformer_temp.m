% Rough Approximation of Transformer Core Temperature
%
% Authors: Julia Estrin
% Date:    03-03-2026
%
%
%   Thermal Model:
%     1. Power dissipation split: 50% in center post, 50% in outer legs/yokes
%     2. Heat path: center post → thermal interface → cold plate → water
%     3. Series thermal resistances: R_ferrite + R_grease + R_plate_water
%     4. Temperature rise: ΔT = P × R_total
%
% Inputs:
%   P_diss        - Total transformer power dissipation (core + copper) [W]
%   A_centerpost  - Center post cross-sectional area [m²]
%   l_centerpost  - Center post length (height) [m]
%   R_plate       - Cold plate thermal resistance (plate to water) [°C/W]
%   A_plate       - Cold plate contact area [m²]
%   T_water       - Water coolant temperature [°C]
%   sig_grease    - Thermal conductivity of interface grease [W/(m·K)]
%   d_grease      - Thermal grease layer thickness [m]
%
% Outputs:
%   T_tx          - Estimated peak transformer temperature (center post) [°C]
%
% Assumptions:
%   - 50/50 power split between center post and outer structure
%   - 1D heat flow from center post to cold plate
%   - Uniform temperature across center post cross-section
%   - Cold plate covers entire center post area
%   - Steady-state thermal conditions
%
% References:
%   Ferrite thermal properties: https://fair-rite.com/tables/
%
% Usage Example:
%   T_tx = calculate_transformer_temp(15, 50e-6, 4e-3, 0.1, 100e-6, 40, 3.5, 100e-6);

function T_tx = calculate_transformer_temp(P_diss, A_centerpost, l_centerpost, R_plate, A_plate, T_water, sig_grease, d_grease)

    %% Power Distribution Assumption
    % Assume half of total power dissipation occurs in the center post,
    % with the remaining half distributed across outer legs and yokes.
    % This is a conservative estimate
    P_centerpost = P_diss / 2;              % [W] Power dissipated in center post

    %% Material Thermal Properties
    % Ferrite thermal conductivity (Fair-Rite typical values for MnZn ferrites)
    % Range: 3.5-4.3 W/(m·K) for typical power ferrites
    % Using conservative value: 3.5 W/(m·K)
    sig_ferrite = 0.1 * 35;                 % [W/(m·K)] Thermal conductivity of ferrite
    rho_ferrite = 1 / sig_ferrite;          % [(m·K)/W] Thermal resistivity of ferrite
    
    % Thermal interface material (grease/paste)
    rho_grease  = 1 / sig_grease;           % [(m·K)/W] Thermal resistivity of grease

    %% Thermal Resistance Network
    % Heat flows vertically through the center post, then through the
    % thermal interface layer, then into the cold plate, and finally to water.
    
    % --- Ferrite Conduction Resistance ---
    % Heat conducts through half the center post height (assuming single-sided cooling from bottom).
    % Thermal resistance: R = ρ·L / A
    R_ferrite = (rho_ferrite * l_centerpost / 2) / A_centerpost;  % [K/W]
    
    % --- Thermal Interface Resistance ---
    % Heat conducts through the thin grease layer between ferrite and cold plate.
    R_grease = (rho_grease * d_grease) / A_centerpost;            % [K/W]
    
    % --- Cold Plate to Water Resistance ---
    % The cold plate has its own thermal resistance R_plate [K/W] defined
    % for a reference area A_plate. Scale this to the actual center post area.
    % If A_centerpost < A_plate, effective resistance increases.
    R_plate_water = R_plate * A_plate / A_centerpost;             % [K/W]
    
    % --- Total Thermal Resistance (Series Path) ---
    R_total = R_ferrite + R_grease + R_plate_water;               % [K/W]

    %% Temperature Calculation
    % Temperature rise above water: ΔT = P × R_total
    % Peak temperature: T = T_water + ΔT
    T_centerpost = P_centerpost * R_total + T_water;              % [°C]
    
    % Return peak temperature (center post is assumed to be the hottest point)
    T_tx = T_centerpost;                                          % [°C]

end