function [results] = optimize_transformer_fmincon(options)
% OPTIMIZE_TRANSFORMER_FMINCON Optimize transformer design using fmincon
%
% Uses gradient-based optimization to minimize losses, much faster than grid sweep
%
% Usage:
%   results = optimize_transformer_fmincon()              % Use defaults
%   results = optimize_transformer_fmincon(options)       % Custom options
%
% Options (optional struct):
%   .a_bounds    - [min, max] bounds for parameter a [m]
%   .b_bounds    - [min, max] bounds for parameter b [m]
%   .w_bounds    - [min, max] bounds for parameter w [m]
%   .a0, .b0, .w0 - Initial guesses [m]
%   .optimize_for - 'Pc' (core loss) or 'Ptot' (total loss)
%   .algorithm   - 'fmincon', 'fminsearch', 'ga', 'particleswarm'
%   .multi_start - Number of random starting points (default: 1)
%   .plot_results - true/false
%
% Returns:
%   results - struct with optimal parameters, losses, and optimization info

    %% Default options
    if nargin < 1
        options = struct();
    end
    
    % Parameter bounds [min, max] in meters
    if ~isfield(options, 'a_bounds'), options.a_bounds = [20e-3, 1000e-3]; end
    if ~isfield(options, 'b_bounds'), options.b_bounds = [5e-3, 1005e-3]; end
    if ~isfield(options, 'w_bounds'), options.w_bounds = [5e-3, 1005e-3]; end
    
    % Initial guess (center of bounds if not specified)
    if ~isfield(options, 'a0'), options.a0 = mean(options.a_bounds); end
    if ~isfield(options, 'b0'), options.b0 = mean(options.b_bounds); end
    if ~isfield(options, 'w0'), options.w0 = mean(options.w_bounds); end
    
    if ~isfield(options, 'optimize_for'), options.optimize_for = 'Ptot'; end
    if ~isfield(options, 'algorithm'), options.algorithm = 'fmincon'; end
    if ~isfield(options, 'multi_start'), options.multi_start = 1; end
    if ~isfield(options, 'plot_results'), options.plot_results = true; end
    if ~isfield(options, 'Bmax_limit'), options.Bmax_limit = 0.35; end

    %% Fixed parameters
    params.material_name = 'F80';
    params.Vo = 48;
    params.f = 1e6;
    params.Np = 4;
    params.L = 3e-6;
    params.u0 = 4*pi*1e-7;
    params.s_ct = 1e-3;
    params.h_pcb = 1.6e-3;
    params.rho_cu = 2.2e-8;
    params.sigma_cu = 1 / params.rho_cu;
    params.t_cu_pri = 2 * 35e-6;
    params.t_cu_sec = 2 * 35e-6;
    params.stackup = '3layer';
    params.I = 20;
    params.Bmax_limit = options.Bmax_limit;

    % Get core material
    material = get_core_material(params.material_name);
    params.k = material.k;
    params.beta = material.beta;
    params.alpha = material.alpha;
    params.uc = material.uc;

    %% Setup optimization
    x0 = [options.a0, options.b0, options.w0];
    lb = [options.a_bounds(1), options.b_bounds(1), options.w_bounds(1)];
    ub = [options.a_bounds(2), options.b_bounds(2), options.w_bounds(2)];

    % Objective function
    obj_fun = @(x) objective(x, params, options.optimize_for);
    
    % Nonlinear constraints (Bmax limit, geometry validity)
    nonlcon = @(x) constraints(x, params);

    %% Run optimization
    fprintf('Starting optimization using %s...\n', options.algorithm);
    fprintf('Optimizing for: %s\n', options.optimize_for);
    fprintf('Bounds: a=[%.1f, %.1f]mm, b=[%.1f, %.1f]mm, w=[%.1f, %.1f]mm\n', ...
        lb(1)*1e3, ub(1)*1e3, lb(2)*1e3, ub(2)*1e3, lb(3)*1e3, ub(3)*1e3);
    
    switch lower(options.algorithm)
        case 'fmincon'
            opts = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'Algorithm', 'sqp', ...
                'MaxFunctionEvaluations', 3000, ...
                'OptimalityTolerance', 1e-8, ...
                'StepTolerance', 1e-10, ...
                'FiniteDifferenceType', 'central');
            
            if options.multi_start > 1
                % MultiStart for global optimization
                problem = createOptimProblem('fmincon', ...
                    'objective', obj_fun, ...
                    'x0', x0, ...
                    'lb', lb, ...
                    'ub', ub, ...
                    'nonlcon', nonlcon, ...
                    'options', opts);
                ms = MultiStart('Display', 'iter');
                [x_opt, fval, exitflag, output, solutions] = run(ms, problem, options.multi_start);
            else
                [x_opt, fval, exitflag, output] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, nonlcon, opts);
            end
            
        case 'fminsearch'
            % Unconstrained - use penalty method
            obj_penalized = @(x) objective_penalized(x, params, options.optimize_for, lb, ub);
            opts = optimset('Display', 'iter', 'MaxFunEvals', 5000, 'TolX', 1e-8);
            [x_opt, fval, exitflag, output] = fminsearch(obj_penalized, x0, opts);
            x_opt = max(min(x_opt, ub), lb);  % Enforce bounds
            
        case 'ga'
            % Genetic algorithm for global optimization
            opts = optimoptions('ga', ...
                'Display', 'iter', ...
                'PopulationSize', 100, ...
                'MaxGenerations', 200, ...
                'FunctionTolerance', 1e-8, ...
                'UseParallel', true);
            [x_opt, fval, exitflag, output] = ga(obj_fun, 3, [], [], [], [], lb, ub, nonlcon, opts);
            
        case 'particleswarm'
            % Particle swarm for global optimization
            opts = optimoptions('particleswarm', ...
                'Display', 'iter', ...
                'SwarmSize', 100, ...
                'MaxIterations', 200, ...
                'FunctionTolerance', 1e-8, ...
                'UseParallel', true);
            % PSO doesn't support nonlinear constraints directly, use penalty
            obj_penalized = @(x) objective_with_constraints(x, params, options.optimize_for);
            [x_opt, fval, exitflag, output] = particleswarm(obj_penalized, 3, lb, ub, opts);
            
        case 'surrogateopt'
            % Surrogate optimization for expensive functions
            opts = optimoptions('surrogateopt', ...
                'Display', 'iter', ...
                'MaxFunctionEvaluations', 500, ...
                'UseParallel', true);
            [x_opt, fval, exitflag, output] = surrogateopt(obj_fun, lb, ub, opts);
            
        otherwise
            error('Unknown algorithm: %s', options.algorithm);
    end

    %% Extract results
    a_opt = x_opt(1);
    b_opt = x_opt(2);
    w_opt = x_opt(3);
    
    % Calculate full loss breakdown at optimal point
    [~, loss_data] = objective(x_opt, params, options.optimize_for);

    %% Store results
    results.a_opt = a_opt;
    results.b_opt = b_opt;
    results.w_opt = w_opt;
    results.Pc_opt = loss_data.Pc;
    results.Ptot_opt = loss_data.Ptot;
    results.P_pri_opt = loss_data.P_pri;
    results.P_sec_opt = loss_data.P_sec;
    results.Bmax_opt = loss_data.Bmax;
    results.l_g1_opt = loss_data.l_g1;
    results.Vc_opt = loss_data.Vc;
    results.Ac_opt = loss_data.Ac;
    results.fval = fval;
    results.exitflag = exitflag;
    results.output = output;
    results.params = params;
    results.options = options;

    %% Display results
    fprintf('\n========== OPTIMIZATION RESULTS ==========\n');
    fprintf('Algorithm: %s\n', options.algorithm);
    fprintf('Exit flag: %d\n', exitflag);
    fprintf('\nOptimal Parameters:\n');
    fprintf('  a = %.3f mm\n', a_opt * 1e3);
    fprintf('  b = %.3f mm\n', b_opt * 1e3);
    fprintf('  w = %.3f mm\n', w_opt * 1e3);
    fprintf('\nLosses at Optimal Point:\n');
    fprintf('  Core Loss (Pc)     = %.4f W\n', loss_data.Pc);
    fprintf('  Primary Loss       = %.4f W\n', loss_data.P_pri);
    fprintf('  Secondary Loss     = %.4f W\n', loss_data.P_sec);
    fprintf('  Total Loss (Ptot)  = %.4f W\n', loss_data.Ptot);
    fprintf('\nDerived Values:\n');
    fprintf('  Bmax = %.4f T (limit: %.2f T)\n', loss_data.Bmax, params.Bmax_limit);
    fprintf('  Air gap (l_g1) = %.4f mm\n', loss_data.l_g1 * 1e3);
    fprintf('  Core Volume (Vc) = %.2f mm³\n', loss_data.Vc * 1e9);
    fprintf('  Core Area (Ac) = %.2f mm²\n', loss_data.Ac * 1e6);
    fprintf('==========================================\n');

    %% Plot results
    if options.plot_results
        plot_sensitivity(results, params, options);
    end
end

%% Objective function
function [f, loss_data] = objective(x, params, optimize_for)
    a = x(1);
    b = x(2);
    w = x(3);
    
    loss_data = calculate_all_losses(a, b, w, params);
    
    if ~loss_data.valid
        f = 1e10;  % Large penalty for invalid designs
        return;
    end
    
    if strcmp(optimize_for, 'Pc')
        f = loss_data.Pc;
    else
        f = loss_data.Ptot;
    end
end

%% Objective with penalty for unconstrained solvers
function f = objective_penalized(x, params, optimize_for, lb, ub)
    % Bound penalty
    penalty = 0;
    for i = 1:3
        if x(i) < lb(i)
            penalty = penalty + 1e6 * (lb(i) - x(i))^2;
        elseif x(i) > ub(i)
            penalty = penalty + 1e6 * (x(i) - ub(i))^2;
        end
    end
    
    [f, ~] = objective(x, params, optimize_for);
    f = f + penalty;
end

%% Objective with constraint penalty for PSO
function f = objective_with_constraints(x, params, optimize_for)
    [f, loss_data] = objective(x, params, optimize_for);
    
    % Add penalty for constraint violations
    if loss_data.valid
        if loss_data.Bmax > params.Bmax_limit
            f = f + 1e6 * (loss_data.Bmax - params.Bmax_limit)^2;
        end
    end
end

%% Nonlinear constraints for fmincon
function [c, ceq] = constraints(x, params)
    a = x(1);
    b = x(2);
    w = x(3);
    
    loss_data = calculate_all_losses(a, b, w, params);
    
    % Inequality constraints: c <= 0
    c = [];
    
    % Bmax constraint
    if loss_data.valid
        c(1) = loss_data.Bmax - params.Bmax_limit;  % Bmax <= Bmax_limit
    else
        c(1) = 1;  % Invalid design
    end
    
    % Minimum winding width
    w_winding = w - 2*params.s_ct;
    c(2) = params.s_ct - w_winding;  % w_winding >= s_ct (minimum viable winding)
    
    % Positive core volume (implicit from geometry)
    c(3) = -loss_data.Vc;  % Vc > 0
    
    % Equality constraints
    ceq = [];
end

%% Calculate all losses
function loss_data = calculate_all_losses(a, b, w, params)
    loss_data = struct();
    loss_data.valid = true;
    
    % Derived geometry
    w_winding = w - 2*params.s_ct;
    
    if w_winding <= 0
        loss_data.valid = false;
        loss_data.Pc = Inf;
        loss_data.Ptot = Inf;
        loss_data.Bmax = Inf;
        loss_data.Vc = 0;
        loss_data.Ac = 0;
        loss_data.l_g1 = 0;
        loss_data.P_pri = 0;
        loss_data.P_sec = 0;
        return;
    end
    
    h_w = w/4 + params.h_pcb;
    w_tot = 4*w + 3*b;
    l = a + 2*w_winding + 2*params.s_ct;
    h = h_w + b;
    Ac = a * b;
    
    l_g1 = params.Np^2 * params.u0 * Ac / (8 * params.L);
    Vc = w_tot * l * h - 2 * w * h_w * a;
    
    if Vc <= 0 || Ac <= 0
        loss_data.valid = false;
        loss_data.Pc = Inf;
        loss_data.Ptot = Inf;
        loss_data.Bmax = Inf;
        loss_data.Vc = 0;
        loss_data.Ac = Ac;
        loss_data.l_g1 = l_g1;
        loss_data.P_pri = 0;
        loss_data.P_sec = 0;
        return;
    end
    
    Bmax = 2 * params.Vo / (4 * params.f * Ac);
    Pc = Vc * params.k * params.f^(params.alpha) * Bmax^(params.beta);
    
    % Winding losses
    Rdc_pri = calculate_Rdc_rectangle(100, ...
                b + 2*params.s_ct, b + 2*params.s_ct + 2*w_winding, ...
                a + 2*params.s_ct, l, ...
                params.sigma_cu, params.t_cu_pri);
    
    Rdc_sec = calculate_Rdc_rectangle(100, ...
                b + 2*params.s_ct, b + 2*params.s_ct + 2*w_winding, ...
                a + 2*params.s_ct, l, ...
                params.sigma_cu, params.t_cu_sec);
    
    [P_pri, P_sec] = calculate_Pcond( ...
            params.I, params.f, ...
            Rdc_pri, Rdc_sec, ...
            params.t_cu_pri, params.t_cu_sec, ...
            params.rho_cu, params.u0, params.stackup);
    
    loss_data.Pc = Pc;
    loss_data.P_pri = P_pri;
    loss_data.P_sec = P_sec;
    loss_data.Ptot = P_pri + P_sec + Pc;
    loss_data.Bmax = Bmax;
    loss_data.l_g1 = l_g1;
    loss_data.Vc = Vc;
    loss_data.Ac = Ac;
    loss_data.w_winding = w_winding;
end

%% Plot sensitivity around optimum
function plot_sensitivity(results, params, options)
    figure('Name', 'Optimization Results', 'Position', [100, 100, 1200, 400]);
    
    a_opt = results.a_opt;
    b_opt = results.b_opt;
    w_opt = results.w_opt;
    
    % Sweep each parameter around optimum
    n_pts = 50;
    
    % a sweep
    a_vec = linspace(options.a_bounds(1), options.a_bounds(2), n_pts);
    Ptot_a = zeros(size(a_vec));
    for i = 1:n_pts
        [~, ld] = objective([a_vec(i), b_opt, w_opt], params, 'Ptot');
        Ptot_a(i) = ld.Ptot;
    end
    
    subplot(1, 3, 1);
    plot(a_vec*1e3, Ptot_a, 'b-', 'LineWidth', 2);
    hold on;
    plot(a_opt*1e3, results.Ptot_opt, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xlabel('a [mm]');
    ylabel('Total Loss [W]');
    title(sprintf('Sensitivity to a (b=%.2f, w=%.2f mm)', b_opt*1e3, w_opt*1e3));
    grid on;
    
    % b sweep
    b_vec = linspace(options.b_bounds(1), options.b_bounds(2), n_pts);
    Ptot_b = zeros(size(b_vec));
    for i = 1:n_pts
        [~, ld] = objective([a_opt, b_vec(i), w_opt], params, 'Ptot');
        Ptot_b(i) = ld.Ptot;
    end
    
    subplot(1, 3, 2);
    plot(b_vec*1e3, Ptot_b, 'b-', 'LineWidth', 2);
    hold on;
    plot(b_opt*1e3, results.Ptot_opt, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xlabel('b [mm]');
    ylabel('Total Loss [W]');
    title(sprintf('Sensitivity to b (a=%.2f, w=%.2f mm)', a_opt*1e3, w_opt*1e3));
    grid on;
    
    % w sweep
    w_vec = linspace(options.w_bounds(1), options.w_bounds(2), n_pts);
    Ptot_w = zeros(size(w_vec));
    for i = 1:n_pts
        [~, ld] = objective([a_opt, b_opt, w_vec(i)], params, 'Ptot');
        Ptot_w(i) = ld.Ptot;
    end
    
    subplot(1, 3, 3);
    plot(w_vec*1e3, Ptot_w, 'b-', 'LineWidth', 2);
    hold on;
    plot(w_opt*1e3, results.Ptot_opt, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xlabel('w [mm]');
    ylabel('Total Loss [W]');
    title(sprintf('Sensitivity to w (a=%.2f, b=%.2f mm)', a_opt*1e3, b_opt*1e3));
    grid on;
    
    sgtitle('Parameter Sensitivity at Optimum', 'FontSize', 14, 'FontWeight', 'bold');
end