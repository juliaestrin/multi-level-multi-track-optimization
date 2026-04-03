function VIRT3Dfigure_ROUNDCenter_V2(result, material, T_tx, idx)
    figure(1)
    % --- Select subplot based on idx ---
    if idx == ~0
        subplot(2, 1, idx);
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

   
  % h_w = w/4 + h_pcb;
  % w_tot = 2*w + 2*b;
  % l = a + 2*w_winding + 2*s_ct;
 
   

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
        x_center = w_tot / 2;
        [X, Y, Z] = cylinder(b/2, 50);
        Z = Z * h;
        X = X + x_center;
        Y = Y + a / 2;
        surf(X, Y, Z, 'FaceColor', core_color, 'EdgeColor', 'none');
        
        % Add top and bottom caps for the cylinder
        % theta = linspace(0, 2*pi, 51);
        % x_cap = r_centerpost * cos(theta);
        % y_cap = r_centerpost * sin(theta) + w_core/2;
        % 
        % fill3(x_cap, y_cap, ones(size(x_cap)) * h_yoke, core_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
        % fill3(x_cap, y_cap, ones(size(x_cap)) * (h_yoke + w_b), core_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
    else
        % Square center post (box)
        drawBox([b/2 + w, 0, 0], [b, a, b/2+h_w-lg], core_color);
    end

    % --- Winding ---
    

    if strcmp(centerpost_shape, 'round')
        % --- Annular (donut-shaped) winding for round centerpost ---
        r_inner = b/2 + s_ct;
        r_outer = r_inner + w_winding;
        y_center = a/2;
        n_theta = 60;
        theta = linspace(0, 2*pi, n_theta + 1);

        z_winding = h/2-h_pcb/2; 
        
        n_radial = 10;
        r_vals = linspace(r_inner, r_outer, n_radial);
        [Theta, R] = meshgrid(theta, r_vals);
        X_ann = R .* cos(Theta) + x_center;
        Y_ann = R .* sin(Theta) + y_center;
        Z_bottom = ones(size(X_ann)) * (z_winding);
        surf(X_ann, Y_ann, Z_bottom, 'FaceColor', copper_pri, 'EdgeColor', 'none');
        
        Z_top_wind = ones(size(X_ann)) * (z_winding + h_pcb);
        surf(X_ann, Y_ann, Z_top_wind, 'FaceColor', copper_pri, 'EdgeColor', 'none');
        
        [X_cyl, Y_cyl, Z_cyl] = cylinder(r_outer, n_theta);
        Z_cyl = Z_cyl * h_pcb + z_winding;
        Y_cyl = Y_cyl + y_center;
        X_cyl = X_cyl + x_center; 
        surf(X_cyl, Y_cyl, Z_cyl, 'FaceColor', copper_pri, 'EdgeColor', 'none');
        
        [X_cyl_in, Y_cyl_in, Z_cyl_in] = cylinder(r_inner, n_theta);
        Z_cyl_in = Z_cyl_in * h_pcb + z_winding;
        Y_cyl_in = Y_cyl_in + y_center;
        X_cyl_in = X_cyl_in + x_center; 
        surf(X_cyl_in, Y_cyl_in, Z_cyl_in, 'FaceColor', copper_pri, 'EdgeColor', 'none');
    else
        % --- Rectangular winding for square centerpost ---
        drawBox([b/2+s_ct, -w_winding - s_ct,  h/2-h_pcb/2], ...
                [w_tot - b - 2 * s_ct, w_winding , h_pcb], copper_pri);
        drawBox([b/2+s_ct, a + s_ct,  h/2-h_pcb/2], ...
                [w_tot - b - 2 * s_ct, w_winding , h_pcb], copper_pri);
        drawBox([b/2+s_ct, -w_winding - s_ct,  h/2-h_pcb/2], ...
                [w_winding, l , h_pcb], copper_pri);
        drawBox([b + w + b/2+s_ct, -w_winding - s_ct,  h/2-h_pcb/2], ...
                [w_winding, l , h_pcb], copper_pri);
    end

    % --- Formatting ---
    xlabel('X [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Z [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    
    title_str = sprintf(['3D EQ-Core Transformer (%s centerpost) - %s\n' ...
                         'P_v = %.0f kW/m³,  Footprint Area = %.1f mm^2,  Total Loss = %.2f W\n' ...
                         'Copper Loss = %.2f W,  Core Loss = %.2f W\n' ...
                         'Temp = %.2f °C'], ...
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
     % 'b          = %.2f mm\n' ...
        % 'w          = %.2f mm\n' ...
        %result.b * 1e3, ...
        %result.w * 1e3, ...
        
    
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