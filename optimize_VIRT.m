% Optimize VIRT Transformer Design over Pv_max and Window Height
%
% Author:  Julia Estrin
% Date:    02-17-2026
%
% Description:
%   Performs a 2D sweep over a range of maximum core loss densities (Pv_max)
%   and winding window heights (w_height) to find the transformer design
%   that minimizes total loss (core + primary copper + secondary copper).
%
%   For each (Pv_max, w_height) pair, calculate_VIRT_design is called and
%   key results are stored in 2D arrays (rows = Pv_max index, columns =
%   w_height index). The global minimum total loss design is identified and
%   its full result struct is returned alongside the sweep arrays.
%
% Inputs:
%   Pv_max_list   - Vector of core loss density values to sweep [W/m^3]
%   w_height_list - Vector of window height values to sweep [m]
%   design_params - Parameter struct (see calculate_VIRT_design for fields)
%
% Outputs:
%   result - Struct containing:
%     Sweep arrays (rows = Pv_max index, cols = w_height index):
%       .P_core_array   - Core loss [W]
%       .P_pri_array    - Primary copper loss [W]
%       .P_sec_array    - Secondary copper loss [W]
%       .P_total_array  - Total loss [W]
%       .V_total_array  - Core volume [m^3]
%       .Rdc_pri_array  - Primary DC resistance [Ohm]
%       .Rdc_sec_array  - Secondary DC resistance [Ohm]
%       .Bmax_array     - Peak flux density [T]
%     Optimal design scalars:
%       .Pv_max_opt     - Optimal core loss density [W/m^3]
%       .w_height_opt   - Optimal window height [m]
%       .P_total_min    - Minimum total loss [W]
%       .P_core_min     - Core loss at optimum [W]
%       .P_copper_min   - Total copper loss at optimum [W]
%       .Bmax_opt       - Peak flux density at optimum [T]
%       .opt_design     - Full result struct from calculate_VIRT_design at optimum
%
% Usage Example:
%   Pv_max_list   = linspace(10e3, 500e3, 50);
%   w_height_list = linspace(5e-3, 50e-3, 50);
%   result = optimize_VIRT(Pv_max_list, w_height_list, design_params);

function result = optimize_VIRT(Pv_max_list, w_height_list, design_params)

    %% Input Validation
    assert(~isempty(Pv_max_list),   'Pv_max_list must not be empty.');
    assert(~isempty(w_height_list), 'w_height_list must not be empty.');

    %% Preallocate Sweep Arrays
    % Rows correspond to Pv_max index, columns to w_height index
    n_Pv = length(Pv_max_list);
    n_wh = length(w_height_list);

    P_core_array  = zeros(n_Pv, n_wh);
    P_pri_array   = zeros(n_Pv, n_wh);
    P_sec_array   = zeros(n_Pv, n_wh);
    P_total_array = zeros(n_Pv, n_wh);
    V_total_array = zeros(n_Pv, n_wh);
    Rdc_pri_array = zeros(n_Pv, n_wh);
    Rdc_sec_array = zeros(n_Pv, n_wh);
    Bmax_array    = zeros(n_Pv, n_wh);

    %% Parameter Sweep
    for i = 1:n_Pv
        for j = 1:n_wh
           % des = calculate_VIRT_design(Pv_max_list(i), w_height_list(j), design_params);
           des = calculate_VIRT_design_02_27_2026(Pv_max_list(i), w_height_list(j), design_params);


            P_core_array(i,j)  = des.P_core;
            P_pri_array(i,j)   = des.P_pri;
            P_sec_array(i,j)   = des.P_sec;
            P_total_array(i,j) = des.P_total;
            V_total_array(i,j) = des.V_total;
            Rdc_pri_array(i,j) = des.Rdc_pri;
            Rdc_sec_array(i,j) = des.Rdc_sec;
            Bmax_array(i,j)    = des.Bmax;
        end
    end

    %% Find Global Optimum
    [P_total_min, linear_idx]    = min(P_total_array(:));
    [idx_Pv_opt, idx_wh_opt]     = ind2sub(size(P_total_array), linear_idx);

    Pv_max_opt    = Pv_max_list(idx_Pv_opt);
    w_height_opt  = w_height_list(idx_wh_opt);
    Bmax_opt      = Bmax_array(idx_Pv_opt, idx_wh_opt);
    P_core_min    = P_core_array(idx_Pv_opt, idx_wh_opt);
    P_copper_min  = P_pri_array(idx_Pv_opt, idx_wh_opt) + P_sec_array(idx_Pv_opt, idx_wh_opt);

    % Re-run design at optimum to get full geometry and loss breakdown
    opt_design = calculate_VIRT_design_02_27_2026(Pv_max_opt, w_height_opt, design_params);

    %% Package Results
    result = struct( ...
        'P_core_array',  P_core_array,  ...
        'P_pri_array',   P_pri_array,   ...
        'P_sec_array',   P_sec_array,   ...
        'P_total_array', P_total_array, ...
        'V_total_array', V_total_array, ...
        'Rdc_pri_array', Rdc_pri_array, ...
        'Rdc_sec_array', Rdc_sec_array, ...
        'Bmax_array',    Bmax_array,    ...
        'Pv_max_opt',    Pv_max_opt,    ...
        'w_height_opt',  w_height_opt,  ...
        'P_total_min',   P_total_min,   ...
        'P_core_min',    P_core_min,    ...
        'P_copper_min',  P_copper_min,  ...
        'Bmax_opt',      Bmax_opt,      ...
        'opt_design',    opt_design     ...
    );

end