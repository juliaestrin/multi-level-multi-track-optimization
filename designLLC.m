% Design LLC Resonant Converter for FCML VIRT Topology
%
% Author:  Julia Estrin
% Date:    02-27-2026
%
% Description:
%   Designs an LLC resonant tank for a Flying Capacitor Multilevel (FCML)
%   VIRT topology. Determines the transformer turns ratio, resonant
%   components (Lr, Cr, Lm), and operating current based on gain
%   requirements across the full load and line regulation range.
%
%   Key features:
%     - Accounts for FCML frequency doubling (f0 = 2*fsw)
%     - Multi-track secondary (nt tracks)
%     - Constrained Qe_max to ensure peak gain occurs at fn <= fn_min
%     - Validates gain curve stays within Mg_min/Mg_max bounds
%
% Inputs:
%   Vin_nom     - Nominal input voltage [V]
%   Vo_nom      - Nominal output voltage [V]
%   Mg_nom      - Nominal LLC gain (typically 1.0)
%   nt          - Number of secondary tracks
%   percentReg  - Line regulation tolerance (e.g., 0.1 for ±10%)
%   fsw         - FCML switching frequency [Hz]
%   f_per       - Frequency range percentage (e.g., 0.25 for ±25%)
%   Pmax        - Maximum output power [W]
%   Pmin        - Minimum output power [W]
%   Ln          - Inductance ratio Lm/Lr
%
% Outputs:
%   result - Struct containing:
%       .N          - Transformer turns ratio (primary:secondary per track)
%       .Cr         - Resonant capacitance [F]
%       .Lr         - Resonant inductance [H]
%       .Lm         - Magnetizing inductance [H]
%       .Qe_max     - Maximum quality factor
%       .Ir_rms     - Primary RMS resonant current at worst case [A]
%       .f0         - Resonant frequency [Hz]
%       .fn_min     - Normalized minimum frequency
%       .fn_max     - Normalized maximum frequency
%       .Mg_min     - Minimum required gain
%       .Mg_max     - Maximum required gain
%       .Re_min     - Minimum equivalent resistance [Ohm]
%       .Re_max     - Maximum equivalent resistance [Ohm]
%
% Usage Example:
%   result = designLLC(1500, 48, 1.0, 2, 0.1, 500e3, 0.25, 6.25e3, 625, 5);

function result = designLLC(Vin_nom, Vo_nom, Mg_nom, nt, percentReg, fsw, f_per, Pmax, Pmin, Ln)

    %% Operating Frequencies
    % FCML frequency doubling: transformer sees 2x switching frequency
    f0 = 2 * fsw;                  % [Hz] Resonant frequency
    
    % Normalized frequency range
    fn_min = 1 - f_per;            % Normalized minimum frequency
    fn_max = 1 + f_per;            % Normalized maximum frequency

    %% Transformer Turns Ratio
    % N is the turns ratio from primary to one secondary track
    % Factor of 4 accounts for FCML voltage division
    N_nom = Mg_nom * (Vin_nom / 4) / Vo_nom;
    N = ceil(N_nom);               % Round up to nearest integer

    %% Gain Range
    % Calculate required gain at voltage extremes
    Vin_max = Vin_nom * (1 + percentReg);
    Vin_min = Vin_nom * (1 - percentReg);
    
    Mg_min = N * Vo_nom / (Vin_max / 4);
    Mg_max = N * Vo_nom / (Vin_min / 4);

    %% Quality Factor
    % Find maximum Qe such that:
    %   1. Mg(fn_min) = Mg_max
    %   2. Peak gain occurs at fn <= fn_min
    %   3. Mg(fn_max) <= Mg_min
    Qe_max = find_Qe_max(Ln, Mg_max, fn_min, Mg_min, fn_max);
    
    if isnan(Qe_max)
        warning('designLLC:NoValidQe', ...
            'No valid Qe found. Design constraints may be infeasible.');
    end

    %% Equivalent Load Resistance
    % Convert output power to equivalent AC resistance
    RL_max = Vo_nom^2 / Pmin;      % [Ohm] Load at minimum power
    RL_min = Vo_nom^2 / Pmax;      % [Ohm] Load at maximum power
    
    % Referred to primary, accounting for multi-track secondary
    % Factor: nt * (8/pi^2) * (N/nt)^2 = (8*N^2)/(pi^2*nt)
    Re_min = nt * (8 * (N/nt)^2) / (pi^2) * RL_min;
    Re_max = nt * (8 * (N/nt)^2) / (pi^2) * RL_max;

    %% Resonant Tank Design
    % Size Cr based on Qe_max and Re_min (worst case)
    Cr = 1 / (2 * pi * Qe_max * f0 * Re_min);   % [F]
    Lr = 1 / ((2 * pi * f0)^2 * Cr);             % [H]
    Lm = Ln * Lr;                                 % [H]

    %% RMS Current Calculation
    % Worst case: maximum power at minimum frequency
    Io_min = Pmin / Vo_nom;        % [A] Minimum output current
    Io_max = Pmax / Vo_nom;        % [A] Maximum output current
    
    % RMS output current (fundamental component)
    Io_max_rms = (pi / (2 * sqrt(2))) * Io_max;
    
    % Primary-side RMS load current per track
    Ioe_max_rms = Io_max_rms / (nt * (N / nt));
    
    % Magnetizing current at minimum frequency
    f0_min = fn_min * f0;
    Voe = (2 * sqrt(2) / pi) * (N / nt) * Vo_nom;  % [V] RMS fundamental output voltage per track
    Im_rms = Voe / (2 * pi * f0_min * Lm);
    
    % Total primary RMS current
    Ir_rms = sqrt(Im_rms^2 + Ioe_max_rms^2);

    %% Package Results
    result = struct( ...
        'N',       N,       ...
        'Cr',      Cr,      ...
        'Lr',      Lr,      ...
        'Lm',      Lm,      ...
        'Qe_max',  Qe_max,  ...
        'Ir_rms',  Ir_rms,  ...
        'f0',      f0,      ...
        'fn_min',  fn_min,  ...
        'fn_max',  fn_max,  ...
        'Mg_min',  Mg_min,  ...
        'Mg_max',  Mg_max,  ...
        'Re_min',  Re_min,  ...
        'Re_max',  Re_max   ...
    );

end

% =========================================================================
function Qe_max = find_Qe_max(Ln, Mg_max, fn_min, Mg_min, fn_max)
% Find Maximum Quality Factor with Gain Constraints
%
% Finds the largest Qe such that:
%   1. Mg(fn_min) = Mg_max
%   2. Peak gain occurs at fn <= fn_min
%   3. Mg(fn_max) <= Mg_min
%
% Inputs:
%   Ln      - Inductance ratio Lm/Lr
%   Mg_max  - Maximum required gain
%   fn_min  - Normalized minimum frequency
%   Mg_min  - Minimum required gain
%   fn_max  - Normalized maximum frequency
%
% Outputs:
%   Qe_max  - Maximum valid quality factor (NaN if no solution exists)

    %% LLC Gain Equation
    % Mg(fn, Qe) = |Ln*fn^2 / ((Ln+1)*fn^2 - 1 + j*(fn^2 - 1)*fn*Qe*Ln)|
    gain_func = @(fn, Qe) abs((Ln * fn.^2) ./ ...
        ((Ln + 1) * fn.^2 - 1 + 1j * (fn.^2 - 1) .* fn * Qe * Ln));
    
    gain_at_fn_min = @(Qe) gain_func(fn_min, Qe);
    error_func = @(Qe) gain_at_fn_min(Qe) - Mg_max;

    %% Coarse Search for Crossing Points
    % Find all Qe values where Mg(fn_min) = Mg_max
    Qe_sweep = logspace(-2, 2, 1000);  % Sweep from 0.01 to 100
    error_sweep = arrayfun(error_func, Qe_sweep);
    
    % Identify sign changes (zero crossings)
    sign_changes = find(diff(sign(error_sweep)) ~= 0);
    
    if isempty(sign_changes)
        Qe_max = NaN;
        return;
    end

    %% Refine Crossing Points
    Qe_candidates = [];
    for idx = sign_changes
        try
            Qe_solution = fzero(error_func, [Qe_sweep(idx), Qe_sweep(idx+1)]);
            if Qe_solution > 0
                Qe_candidates = [Qe_candidates, Qe_solution];
            end
        catch
            continue;
        end
    end
    
    if isempty(Qe_candidates)
        Qe_max = NaN;
        return;
    end

    %% Validate Constraints for Each Candidate
    valid_Qe = [];
    
    for Qe = Qe_candidates
        % --- Constraint 1: Peak gain must occur at fn <= fn_min ---
        % Find peak frequency for this Qe
        fn_search = linspace(0.3, 1.5, 1000);
        Mg_search = arrayfun(@(f) gain_func(f, Qe), fn_search);
        
        [~, max_idx] = max(Mg_search);
        fn_peak = fn_search(max_idx);
        
        % Refine peak location using fminbnd
        try
            if max_idx > 1 && max_idx < length(fn_search)
                fn_left  = fn_search(max(1, max_idx - 10));
                fn_right = fn_search(min(length(fn_search), max_idx + 10));
                [fn_peak, ~] = fminbnd(@(fn) -gain_func(fn, Qe), fn_left, fn_right);
            end
        catch
            % Use coarse estimate if refinement fails
        end
        
        % --- Constraint 2: Mg(fn_max) <= Mg_min ---
        Mg_at_fn_max = gain_func(fn_max, Qe);
        
        % Check both constraints
        if (fn_peak <= fn_min) && (Mg_at_fn_max <= Mg_min)
            valid_Qe = [valid_Qe, Qe];
        end
    end

    %% Return Maximum Valid Qe
    if ~isempty(valid_Qe)
        Qe_max = max(valid_Qe);
    else
        Qe_max = NaN;  % No valid solution found
    end

end