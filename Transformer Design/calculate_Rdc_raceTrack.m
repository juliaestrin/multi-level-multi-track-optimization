% Calculate Equivalent DC Resistance of a Racetrack Planar Winding
%
% Author:  Julia Estrin
% Date:    06-14-2026
%
% Description:
%   Computes the equivalent DC resistance of a racetrack (stadium) planar
%   winding by modeling the total winding as (n+1) concentric racetrack
%   turns, each with linearly increasing dimensions. The per-turn resistances
%   are then combined in parallel to yield the equivalent DC resistance.
%
%   A racetrack turn consists of two straight sections and two semicircular
%   ends of diameter b. The straight-section contribution to the perimeter
%   is (a - b) per side (the chord length excluding the semicircle caps),
%   giving a per-turn perimeter of:
%       L = 2*(a(i) - b(i)) + pi*b(i)
%
%   The wall thickness w is uniform everywhere. The innermost turn
%   centerline sits at (a + dw, b + dw) and the outermost at
%   (a + 2w - dw, b + 2w - dw), where dw = w / (n+1).
%   Both a and b grow by 2*dw per turn, preserving uniform wall width.
%
% Inputs:
%   n       - Index of the outermost turn (loop runs from turn 0 to turn n)
%   a       - Inner straight-side length (innermost turn boundary) [m]
%   b       - Inner end diameter (innermost turn boundary) [m]
%   w       - Uniform wall (winding) thickness [m]
%   sig_cu  - Electrical conductivity of copper [S/m]
%   t_cu    - Copper trace thickness [m]
%
% Outputs:
%   Rdc     - Equivalent DC resistance of the winding [Ohm]
%
% Geometry Notes:
%   - Innermost boundary:  a(1)   = a,       b(1)   = b
%   - Outermost boundary:  a(end) = a + 2*w, b(end) = b + 2*w
%   - Conductor radial width: dw = w / (n + 1)
%   - Turn centerline i:   a(i) = a + (2*i-1)*dw,  b(i) = b + (2*i-1)*dw
%   - Perimeter per turn:  L(i) = 2*(a(i) - b(i)) + pi*b(i)
%   - Cross-section area:  A = dw * t_cu
%   - Per-turn resistance: R(i) = L(i) / (sig_cu * dw * t_cu)
%   - All turns combined in parallel for equivalent resistance
%
% Usage Example:
%   % 10 turns, inner a=40mm, inner b=20mm, wall w=5mm
%   % 2oz copper (70um), sigma = 4.5e7 S/m
%   Rdc = calculate_Rdc_racetrack(10, 40e-3, 20e-3, 5e-3, 4.5e7, 70e-6);

function Rdc = calculate_Rdc_raceTrack(n, a, b, w, sig_cu, t_cu)

% --- Geometry: conductor radial width per turn ---
% dw is the radial width of each conductor strip.
% w is distributed evenly across (n+1) conductors.
dw = w / (n + 1);

% --- Preallocate arrays for performance ---
a_i = zeros(1, n + 1);   % straight-side length centerline per turn [m]
b_i = zeros(1, n + 1);   % end diameter centerline per turn [m]
R   = zeros(1, n + 1);   % resistance per turn [Ohm]

% --- Compute per-turn dimensions and resistances ---
% Iterates over n+1 turns, indexed 1 through n+1 (i = 1 ... n+1).
% Turn centerline dimensions are linearly spaced, with the first turn
% centerline at (a + dw, b + dw) and each subsequent turn stepping
% outward by 2*dw (one dw on each side).
%
% Resistance formula:
%   R = L / (sigma * A)
%   L = 2*(a_i - b_i) + pi*b_i   (racetrack perimeter)
%   A = dw * t_cu                 (cross-section, uniform around perimeter)
for i = 1 : n + 1
    % Centerline dimensions of turn i
    a_i(i) = a + (2*i - 1) * dw;
    b_i(i) = b + (2*i - 1) * dw;

    % Racetrack perimeter and per-turn resistance
    L    = 2 * (a_i(i) - b_i(i)) + pi * b_i(i);
    R(i) = L / (sig_cu * dw * t_cu);
end

% --- Combine resistances ---
% All turns are combined in parallel for the equivalent winding resistance.
Rdc = 1 / sum(1 ./ R);

end