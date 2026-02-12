% Calculate Equivalent DC Resistance of a Rectangular Planar Winding
%
% Author:  Julia Estrin
% Date:    02-12-2026
%
% Description:
%   Computes the equivalent DC resistance of a rectangular planar winding
%   by modeling the total winding as (n+1) series-connected turns, each
%   with linearly increasing side lengths from the innermost turn (a_0, b_0)
%   to the outermost turn (a_n, b_n). The per-turn resistances are then
%   combined in parallel to yield the equivalent DC resistance.
%
%   The loop iterates over (n+1) turns, indexed 0 through n, where n is
%   the index of the outermost turn.
%
% Inputs:
%   n       - Index of the outermost turn (loop runs from turn 0 to turn n)
%   a_0     - Inner winding length in the x-direction [m]
%   a_n     - Outer winding length in the x-direction [m]
%   b_0     - Inner winding length in the y-direction [m]
%   b_n     - Outer winding length in the y-direction [m]
%   sig_cu  - Electrical conductivity of copper [S/m]
%   t_cu    - Copper trace thickness [m]
%
% Outputs:
%   Rdc     - Equivalent DC resistance of the winding [Ohm]
%
% Usage Example:
%   Rdc = calculate_Rdc_summation(5, 0.01, 0.05, 0.01, 0.05, 5.8e7, 35e-6);

function Rdc = calculate_Rdc_summation(n, a_0, a_n, b_0, b_n, sig_cu, t_cu)

    % --- Geometry: conductor width increments per turn ---
    % dx and dy represent the width of each conductor strip in x and y,
    % computed by distributing the total growth evenly across (n+1) turns.
    % NOTE: the divisor uses (n+1) to avoid division by zero when n = 0,
    % and reflects the spacing between (n+2) boundaries for (n+1) strips.
    dx = (a_n - a_0) / (2 * (n + 1));
    dy = (b_n - b_0) / (2 * (n + 1));

    % --- Preallocate arrays for performance ---
    a  = zeros(1, n + 1);   % x-direction half-lengths per turn [m]
    b  = zeros(1, n + 1);   % y-direction half-lengths per turn [m]
    Ra = zeros(1, n + 1);   % x-side resistance per turn [Ohm]
    Rb = zeros(1, n + 1);   % y-side resistance per turn [Ohm]

    % --- Compute per-turn dimensions and resistances ---
    % Iterates over n+1 turns, indexed 0 through n (j = i - 1).
    % Each turn is a rectangle with side lengths increasing linearly.
    % Ra accounts for the two x-direction sides; Rb for the two y-direction sides.
    % Conductor cross-section for x-sides: dy * t_cu
    % Conductor cross-section for y-sides: dx * t_cu
    for i = 1 : n + 1
        j    = i - 1;
        a(i) = a_0 + j * (a_n - a_0) / n;
        b(i) = b_0 + j * (b_n - b_0) / n;
        Ra(i) = 2 * a(i) / (sig_cu * dy * t_cu);
        Rb(i) = 2 * b(i) / (sig_cu * dx * t_cu);
    end

    % --- Combine resistances ---
    % Each turn's total resistance is the sum of its x and y side resistances.
    % All turns are combined in parallel for the equivalent winding resistance.
    Rc      = Ra + Rb;
    Rdc     = 1 / sum(1 ./ Rc);

end