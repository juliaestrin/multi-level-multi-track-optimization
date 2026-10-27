function R = calculate_Rdc_round(r1, r2, sigma, t)

%CALCULATE_RDC_ROUND  DC resistance of a cut annulus, circumferential current.
%
%   R = calculate_Rdc_round(r1, r2, sigma, t)
%
%   Inputs:
%     r1    - inner radius [m]
%     r2    - outer radius [m]   (0 < r1 < r2)
%     sigma - conductivity [S/m]
%     t     - thickness [m]
%   Output:
%     R     - resistance [Ohm]
%
%   Result:
%     R = 2*pi / (sigma * t * log(r2/r1))
%
%   ------------------------------------------------------------------
%   DERIVATION
%   ------------------------------------------------------------------
%   Current flows around the ring (circumferentially), entering and
%   leaving at the two faces of a radial cut.
%
%   1) Slice the annulus into thin rings of radius r and width dr.
%
%   2) Resistance of one ring, using R = L/(sigma*A):
%        Path length:    L = 2*pi*(r + dr/2) = 2*pi*r + pi*dr
%        Cross-section:  A = t*dr
%
%      The pi*dr term in L is a small correction to the finite 2*pi*r
%      (relative size dr/(2r) -> 0), so L -> 2*pi*r. In A, dr is the
%      whole quantity and must be kept.
%
%        dR = 2*pi*r / (sigma*t*dr)
%
%   3) All rings span the same cut, so they share the same voltage
%      and are in PARALLEL. Parallel elements add as conductances:
%
%        dG = 1/dR = (sigma*t / (2*pi)) * dr/r
%
%   4) Integrate from r1 to r2:
%
%        G = (sigma*t/(2*pi)) * integral_{r1}^{r2} (1/r) dr
%          = (sigma*t/(2*pi)) * [ln(r2) - ln(r1)]
%          = (sigma*t/(2*pi)) * ln(r2/r1)
%
%   5) Invert:
%
%        R = 1/G = 2*pi / (sigma * t * ln(r2/r1))
%
%   Note on dropping the dr/2 term:
%     Keeping it gives dG ~ (sigma*t/(2*pi))*(dr/r) - (sigma*t/(4*pi))*(dr^2/r^2).
%     Summing N = (r2-r1)/dr slices, the dr^2 terms total ~dr -> 0,
%     so only the first-order term survives.
%
%   Sanity check (thin ring, r2 - r1 << r1):
%     ln(r2/r1) ~ (r2-r1)/r1, so R ~ 2*pi*r1 / (sigma*t*(r2-r1)),
%     i.e. a straight strip of length 2*pi*r1, width (r2-r1), thickness t.
%   ------------------------------------------------------------------

    if ~(r1 > 0 && r2 > r1)
        error('calculate_Rdc_round:badRadii', 'Require 0 < r1 < r2.');
    end
    R = 2*pi ./ (sigma .* t .* log(r2 ./ r1));
end