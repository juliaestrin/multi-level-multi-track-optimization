function Llk = calc_Llk(design_params, results, h_prepreg, h_core)
% calc_Llk
% Calculates primary-referred leakage inductance for an 8:1 planar transformer
% using the magnetic energy method with radial integration.
%
% Assumptions:
%   - 8 primary layers, each with 1 turn
%   - 2 secondary layers in parallel
%   - each secondary layer carries 4*Ip
%   - leakage field is circumferential: H(r,x) = n(x)*Ip/(2*pi*r)
%   - Llk is referred to the primary
%
% Output:
%   Llk [H]

%% Geometry and material parameters

mu0 = design_params.u0;

t_pri    = design_params.t_cu_pri;
t_sec    = design_params.t_cu_sec;
t_shield = t_pri;   % change this if shield copper thickness differs

r_i = results.b/2 + design_params.s_ct;
r_o = r_i + results.w_winding;

%% Stackup thicknesses

h1  = h_prepreg;
h2  = h_core;
h3  = h_prepreg;
h4  = h_core;
h5  = h_prepreg;
h6  = h_core;
h7  = h_prepreg;
h8  = h_core;
h9  = h_prepreg;
h10 = h_core;
h11 = h_prepreg;

%% Helper functions

% Constant-field insulation/core contribution:
% integral n^2 dx = n^2*h
gap_term = @(n,h) n.^2 .* h;

% Linear-ramp copper contribution:
% integral n(x)^2 dx = t*(a^2 + a*b + b^2)/3
cu_term = @(a,b,t) t .* (a.^2 + a.*b + b.^2) ./ 3;

%% Vertical stack integral: int n(x)^2 dx
%
% Stack assumed from top to bottom:
%
%   Sec top:      0 -> 4
%   h1:          4
%   Shield top:  4 -> 4
%   h2:          4
%   Pri 1:       4 -> 3
%   h3:          3
%   Pri 2:       3 -> 2
%   h4:          2
%   Pri 3:       2 -> 1
%   h5:          1
%   Pri 4:       1 -> 0
%   h6:          0
%   Pri 5:       0 -> 1
%   h7:          1
%   Pri 6:       1 -> 2
%   h8:          2
%   Pri 7:       2 -> 3
%   h9:          3
%   Pri 8:       3 -> 4
%   h10:         4
%   Shield bot:  4 -> 4
%   h11:         4
%   Sec bot:     4 -> 0

stack_integral = ...
    cu_term(0,4,t_sec)       + ...
    gap_term(4,h1)           + ...
    cu_term(4,4,t_shield)    + ...
    gap_term(4,h2)           + ...
    cu_term(4,3,t_pri)       + ...
    gap_term(3,h3)           + ...
    cu_term(3,2,t_pri)       + ...
    gap_term(2,h4)           + ...
    cu_term(2,1,t_pri)       + ...
    gap_term(1,h5)           + ...
    cu_term(1,0,t_pri)       + ...
    gap_term(0,h6)           + ...
    cu_term(0,-1,t_pri)       + ...
    gap_term(-1,h7)           + ...
    cu_term(-1,-2,t_pri)       + ...
    gap_term(-2,h8)           + ...
    cu_term(-2,-3,t_pri)       + ...
    gap_term(-3,h9)           + ...
    cu_term(-3,-4,t_pri)       + ...
    gap_term(-4,h10)          + ...
    cu_term(-4,-4,t_shield)    + ...
    gap_term(-4,h11)          + ...
    cu_term(-4,0,t_sec);

%% Radial integration factor

radial_factor = 2*pi*mu0 / log(r_o/r_i);

%% Primary-referred leakage inductance

Llk = radial_factor * stack_integral;

end