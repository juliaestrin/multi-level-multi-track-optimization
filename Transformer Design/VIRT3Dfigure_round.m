function VIRT3Dfigure_round(result, material, T_tx, idx)
    figure(1)
    % --- Select subplot based on idx ---
    if idx == 0
    elseif idx == 1
        subplot(2, 2, 1);
    elseif idx == 2
        subplot(2, 2, 2);
    elseif idx == 3
        subplot(2, 2, 3);
    elseif idx == 4
        subplot(2, 2, 4);
    end
    cla;
    hold on;
    axis equal;
    view(35, 15);
    grid on;
    grid minor;
    box on;

    % --- Colors ---
    core_color = [0.5, 0.5, 0.8];
    copper_pri = [0.85, 0.45, 0.01];

    a          = result.a          * 1e3;
    b          = result.b          * 1e3;
    w          = result.w          * 1e3;
    Ac         = result.Ac         * 1e6;
    s_ct       = result.s_ct       * 1e3;
    h_pcb      = result.h_pcb      * 1e3;
    lg         = result.lg         * 1e3;
    w_winding  = result.w_winding  * 1e3;
    h_w        = result.h_w        * 1e3;
    w_tot      = result.w_tot      * 1e3;
    h          = result.h          * 1e3;

    centerpost_shape = result.centerpost_shape;
    gap_loc          = result.gap_loc;  % 'center' or 'all'

    w_leg    = 1/2 * Ac / a;
    x_center = w_tot / 2;

    % --- Gap sits just below the top yoke ---
    % Legs run from z=0 to z=(w_leg + h_w), then gap, then top yoke.
    % leg_top_z: top of the gapped leg core piece (bottom of gap)
    % gap sits between leg_top_z and the top yoke (w_leg + h_w)
    leg_top_z = w_leg + h_w - lg;   % bottom of gap = top of lower leg piece

    % --- Bottom Yoke ---
    drawBox([0, 0, 0], [w_tot, a, w_leg], core_color);

    % --- Top Yoke ---
    drawBox([0, 0, w_leg + h_w], [w_tot, a, w_leg], core_color);

    % --- Left Leg ---
    if strcmp(gap_loc, 'all')
        % Lower piece (full height minus gap)
        drawBox([0, 0, w_leg], [w_leg, a, h_w - lg], core_color);
        % gap is left empty between leg_top_z and (w_leg + h_w)
    else
        % No gap on outer legs — draw full height
        drawBox([0, 0, w_leg], [w_leg, a, h_w], core_color);
    end

    % --- Right Leg ---
    if strcmp(gap_loc, 'all')
        drawBox([w_tot - w_leg, 0, w_leg], [w_leg, a, h_w - lg], core_color);
    else
        drawBox([w_tot - w_leg, 0, w_leg], [w_leg, a, h_w], core_color);
    end

    % --- Center Post ---
    % Always has a gap at the top; draw cylinder from z=0 to leg_top_z
    if strcmp(centerpost_shape, 'round')
        % Lower portion of center post (below gap)
        [X, Y, Z] = cylinder(b/2, 50);
        Z = Z * leg_top_z;
        X = X + x_center;
        Y = Y + a / 2;
        surf(X, Y, Z, 'FaceColor', core_color, 'EdgeColor', 'none');
    else
        drawBox([x_center - b/2, 0, 0], [b, a, leg_top_z], core_color);
    end

    % --- Winding (annular donut around center post) ---
    r_inner  = b/2 + s_ct;
    r_outer  = r_inner + w_winding;
    y_center = a / 2;
    n_theta  = 60;
    theta    = linspace(0, 2*pi, n_theta + 1);

    z_bottom = h/2 - h_pcb/2;
    z_top    = h/2 + h_pcb/2;

    n_radial = 10;
    r_vals   = linspace(r_inner, r_outer, n_radial);
    [Theta, R] = meshgrid(theta, r_vals);
    X_ann = R .* cos(Theta) + x_center;
    Y_ann = R .* sin(Theta) + y_center;

    % Bottom face
    surf(X_ann, Y_ann, ones(size(X_ann)) * z_bottom, ...
        'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);
    % Top face
    surf(X_ann, Y_ann, ones(size(X_ann)) * z_top, ...
        'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);

    % Outer cylindrical wall
    [X_cyl, Y_cyl, Z_cyl] = cylinder(r_outer, n_theta);
    Z_cyl = Z_cyl * h_pcb + z_bottom;
    surf(X_cyl + x_center, Y_cyl + y_center, Z_cyl, ...
        'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);

    % Inner cylindrical wall
    [X_cyl_in, Y_cyl_in, Z_cyl_in] = cylinder(r_inner, n_theta);
    Z_cyl_in = Z_cyl_in * h_pcb + z_bottom;
    surf(X_cyl_in + x_center, Y_cyl_in + y_center, Z_cyl_in, ...
        'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);

    % --- Formatting ---
    xlabel('X [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Z [mm]', 'FontSize', 12, 'FontWeight', 'bold');

    title_str = sprintf(['3D EQ-Core Transformer - %s (%s centerpost) - %s\n' ...
                         'P_v = %.0f kW/m³,  Footprint Area = %.1f mm^2,  Total Loss = %.2f W\n' ...
                         'Copper Loss = %.2f W,  Core Loss = %.2f W\n' ...
                         'Temp = %.2f °C,  Gap = %.3f mm (%s legs)'], ...
                        result.topology, ...
                        centerpost_shape, ...
                        material.name, ...
                        result.Pv, ...
                        result.A_footprint * 1e6, ...
                        result.P_total, ...
                        result.P_pri + result.P_sec, ...
                        result.P_core, ...
                        T_tx, ...
                        result.lg * 1e3, ...
                        gap_loc);
    title(title_str, 'FontSize', 13, 'FontWeight', 'bold');

    % --- Info textbox ---
    info_str = sprintf([ ...
        '-- Losses --\n' ...
        'P_pri    = %.3f W\n' ...
        'P_sec    = %.3f W\n' ...
        'P_core   = %.3f W\n' ...
        'P_total  = %.3f W\n' ...
        '\n-- Geometry --\n' ...
        'w_tot [X]    = %.2f mm\n' ...
        'w_winding    = %.2f mm\n' ...
        'l_tot [Y]    = %.2f mm\n' ...
        'l_core       = %.2f mm\n' ...
        'h_core [Z]   = %.2f mm\n' ...
        'l_gap        = %.3f mm\n' ...
        'gap loc      = %s'], ...
        result.P_pri, ...
        result.P_sec, ...
        result.P_core, ...
        result.P_total, ...
        result.w_tot   * 1e3, ...
        result.w_winding * 1e3, ...
        result.l       * 1e3, ...
        result.a       * 1e3, ...
        result.h       * 1e3, ...
        result.lg      * 1e3, ...
        gap_loc);

    % Place info box in a fixed-size panel anchored to the figure,
    % bypassing the annotation clipping issue entirely.
    % Units are normalized [left, bottom, width, height].
    info_panel = uipanel('Units', 'normalized', ...
        'Position',        [0.74, 0.25, 0.25, 0.45], ...
        'BackgroundColor', 'w', ...
        'BorderType',      'line', ...
        'HighlightColor',  'k', ...
        'FontSize',        10);
    uicontrol(info_panel, ...
        'Style',           'text', ...
        'Units',           'normalized', ...
        'Position',        [0.03, 0.02, 0.94, 0.96], ...
        'String',          info_str, ...
        'FontSize',        10, ...
        'FontName',        'FixedWidth', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w');

    rotate3d on;
    hold off;
end