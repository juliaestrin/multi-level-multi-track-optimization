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
    cla;  % Clear current axes instead of clf
    hold on;
    axis equal;
    view(35, 15);
    grid on;
    grid minor;
    box on;

    % --- Colors ---
    core_color = [0.5, 0.5, 0.8];
    copper_pri = [0.85, 0.45, 0.01];

   a = result.a * 1e3; 
   b = result.b * 1e3; 
   w = result.w * 1e3; 
   Ac = result.Ac * 1e6;
   s_ct = result.s_ct * 1e3; 
   h_pcb = result.h_pcb * 1e3; 
   lg = result.lg * 1e3; 
   w_winding = result.w_winding*1e3;
   h_w = result.h_w * 1e3; 
   w_tot = result.w_tot * 1e3; 
   h = result.h * 1e3; 
   
   centerpost_shape = result.centerpost_shape;

   w_leg = 1/2 * Ac / a; 
   x_center = w_tot / 2;

    % --- Bottom Yoke ---
    drawBox([0, 0, 0], [w_tot, a, w_leg], core_color);

    % --- Top Yoke ---
    drawBox([0, 0, w_leg + h_w], [w_tot, a, w_leg], core_color);

    % --- Left Leg ---
    drawBox([0, 0, 0], [w_leg, a, h], core_color);

    % --- Right Leg ---
    drawBox([w_tot - w_leg, 0, 0], [w_leg, a, h], core_color);

    % --- Center Post (shape-dependent) ---
    if strcmp(centerpost_shape, 'round')
        % Cylindrical center post
        [X, Y, Z] = cylinder(b/2, 50);
        Z = Z * h;
        X = X + x_center;
        Y = Y + a / 2;
        surf(X, Y, Z, 'FaceColor', core_color, 'EdgeColor', 'none');
    else
        % Square center post (box)
        drawBox([b/2 + w, 0, 0], [b, a, b/2+h_w-lg], core_color);
    end

    % --- Annular (donut-shaped) winding for round centerpost ---
    r_inner = b/2 + s_ct;
    r_outer = r_inner + w_winding;
    y_center = a/2;
    n_theta = 60;
    theta = linspace(0, 2*pi, n_theta + 1);

    z_bottom = h/2 - h_pcb/2;  % Bottom of winding
    z_top = h/2 + h_pcb/2;     % Top of winding
    
    % --- Top and bottom annular faces ---
    n_radial = 10;
    r_vals = linspace(r_inner, r_outer, n_radial);
    [Theta, R] = meshgrid(theta, r_vals);
    X_ann = R .* cos(Theta) + x_center;
    Y_ann = R .* sin(Theta) + y_center;
    
    % Bottom face
    Z_bottom_face = ones(size(X_ann)) * z_bottom;
    surf(X_ann, Y_ann, Z_bottom_face, 'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);
    
    % Top face
    Z_top_face = ones(size(X_ann)) * z_top;
    surf(X_ann, Y_ann, Z_top_face, 'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);
    
    % --- Outer cylindrical wall ---
    [X_cyl, Y_cyl, Z_cyl] = cylinder(r_outer, n_theta);
    Z_cyl = Z_cyl * h_pcb + z_bottom;  % Scale and shift Z
    X_cyl = X_cyl + x_center;
    Y_cyl = Y_cyl + y_center;
    surf(X_cyl, Y_cyl, Z_cyl, 'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);
    
    % --- Inner cylindrical wall ---
    [X_cyl_in, Y_cyl_in, Z_cyl_in] = cylinder(r_inner, n_theta);
    Z_cyl_in = Z_cyl_in * h_pcb + z_bottom;  % Scale and shift Z
    X_cyl_in = X_cyl_in + x_center;
    Y_cyl_in = Y_cyl_in + y_center;
    surf(X_cyl_in, Y_cyl_in, Z_cyl_in, 'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);

    % --- Formatting ---
    xlabel('X [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Z [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    
    title_str = sprintf(['3D EQ-Core Transformer - %s (%s centerpost) - %s\n' ...
                         'P_v = %.0f kW/m³,  Footprint Area = %.1f mm^2,  Total Loss = %.2f W\n' ...
                         'Copper Loss = %.2f W,  Core Loss = %.2f W\n' ...
                         'Temp = %.2f °C'], ...
                        result.topology, ...
                        centerpost_shape, ...
                        material.name, ...
                        result.Pv, ...
                        result.A_footprint * 1e6, ...
                        result.P_total, ...
                        result.P_pri + result.P_sec, ...
                        result.P_core, ...
                        T_tx);
    title(title_str, 'FontSize', 13, 'FontWeight', 'bold');

    % --- Add textbox with design parameters ---
    info_str = sprintf([ ...
        '\\bf{Losses}\\rm\n' ...
        'P_{pri}    = %.3f W\n' ...
        'P_{sec}    = %.3f W\n' ...
        'P_{core}   = %.3f W\n' ...
        'P_{total}  = %.3f W\n' ...
        '\n\\bf{Geometry}\\rm\n' ... 
        'w_{tot} [X]       = %.2f mm\n' ...
        'w_{winding}       = %.2f mm\n' ...
        'l_{tot}[Y]        = %.2f mm\n' ...
        'l_{core}          = %.2f mm\n' ...
        'h_{core} [Z]      = %.2f mm'], ...
        result.P_pri, ...
        result.P_sec, ...
        result.P_core, ...
        result.P_total, ...
        result.w_tot * 1e3, ...
        result.w_winding * 1e3, ...
        result.l * 1e3, ...
        result.a * 1e3, ...
        result.h * 1e3);
    
    annotation('textbox', [0.75, 0.3, 0.2, 0.4], ...
        'String', info_str, ...
        'FontSize', 10, ...
        'FontName', 'FixedWidth', ...
        'EdgeColor', 'k', ...
        'BackgroundColor', 'w', ...
        'FitBoxToText', 'on', ...
        'Interpreter', 'tex');

    rotate3d on;
    hold off;

end