% Calculate Equivalent DC Resistance of a Round (Annular) Planar Winding
%
% Author:  Julia Estrin
% Date:    03-27-2026
%
% Description:
%   Computes the equivalent DC resistance of a round (annular) planar winding
%   by modeling the total winding as (n+1) series-connected concentric circular
%   turns, each with linearly increasing radius from the innermost turn (r_0)
%   to the outermost turn (r_n). The per-turn resistances are then combined
%   in parallel to yield the equivalent DC resistance.
%
%   Each turn is modeled as a circular loop with circumference 2*pi*r.
%   The conductor cross-section is (radial width dr) x (thickness t_cu).
%
%   The loop iterates over (n+1) turns, indexed 0 through n, where n is
%   the index of the outermost turn.
%
% Inputs:
%   n       - Index of the outermost turn (loop runs from turn 0 to turn n)
%   r_0     - Inner winding radius [m]
%   r_n     - Outer winding radius [m]
%   sig_cu  - Electrical conductivity of copper [S/m]
%   t_cu    - Copper trace thickness [m]
%
% Outputs:
%   Rdc     - Equivalent DC resistance of the winding [Ohm]
%
% Geometry Notes:
%   - Each turn is a circular trace at radius r(i)
%   - Trace length per turn: L(i) = 2 * pi * r(i)
%   - Trace cross-section area: A = dr * t_cu
%   - Per-turn resistance: R(i) = L(i) / (sig_cu * A) = 2*pi*r(i) / (sig_cu * dr * t_cu)
%   - All turns combined in parallel for equivalent resistance
%
% Usage Example:
%   % 100 turns from r=5mm to r=25mm, 2oz copper (70um), at 100°C
%   Rdc = calculate_Rdc_round(100, 5e-3, 25e-3, 4.5e7, 70e-6);

function Rdc = calculate_Rdc_round(n, r_0, r_n, sig_cu, t_cu)

    % --- Input validation ---
    if n < 0
        error('n must be non-negative');
    end
    if r_n <= r_0
        error('Outer radius r_n must be greater than inner radius r_0');
    end
    if r_0 <= 0
        error('Inner radius r_0 must be positive');
    end

    % --- Geometry: conductor radial width per turn ---
    % dr represents the radial width of each conductor strip,
    % computed by distributing the total radial growth evenly across (n+1) turns.
    % Total radial span is (r_n - r_0), divided among (n+1) conductors.
    dr = (r_n - r_0) / (n + 1);

    % --- Preallocate arrays for performance ---
    r = zeros(1, n + 1);    % radius of each turn centerline [m]
    R = zeros(1, n + 1);    % resistance per turn [Ohm]

    % --- Compute per-turn dimensions and resistances ---
    % Iterates over n+1 turns, indexed 0 through n (j = i - 1).
    % Each turn is a circular trace at radius r(i).
    % Turn centerline radii are spaced from (r_0 + dr/2) to (r_n - dr/2).
    %
    % Resistance formula:
    %   R = L / (sigma * A)
    %   L = 2 * pi * r       (circumference)
    %   A = dr * t_cu        (cross-section area)
    %   R = 2 * pi * r / (sig_cu * dr * t_cu)

    for i = 1 : n + 1
        j = i - 1;
        
        % Centerline radius of turn j (linearly spaced)
        % First turn at r_0 + dr/2, last turn at r_n - dr/2
        r(i) = r_0 + dr/2 + j * dr;
        
        % Resistance of this circular turn
        R(i) = 2 * pi * r(i) / (sig_cu * dr * t_cu);
    end

    % --- Combine resistances ---
    % All turns are combined in parallel for the equivalent winding resistance.
    Rdc = 1 / sum(1 ./ R);

end