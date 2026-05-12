function Lr_3DFigure(result)
    figure(2)
    
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
    s_ct       = result.s_ct       * 1e3;
    h_pcb      = result.h_pcb      * 1e3;
    lg         = result.lg         * 1e3;
    w_winding  = result.w_winding  * 1e3;
    h_w        = result.h_w        * 1e3;
    w_tot      = result.w_tot      * 1e3;
    h          = result.h          * 1e3;

    shape    = result.shape;
    gap_loc  = result.gap_loc;

    w_leg    = result.w_leg * 1e3;
    x_center = w_tot / 2;
    y_center = a / 2;

    leg_top_z = w_leg + h_w - lg;

    % --- Bottom Yoke ---
    drawBox([0, 0, 0], [w_tot, a, w_leg], core_color);

    % --- Top Yoke ---
    drawBox([0, 0, w_leg + h_w], [w_tot, a, w_leg], core_color);

    % --- Left Leg ---
    if strcmp(gap_loc, 'all')
        drawBox([0, 0, w_leg], [w_leg, a, h_w - lg], core_color);
    else
        drawBox([0, 0, w_leg], [w_leg, a, h_w], core_color);
    end

    % --- Right Leg ---
    if strcmp(gap_loc, 'all')
        drawBox([w_tot - w_leg, 0, w_leg], [w_leg, a, h_w - lg], core_color);
    else
        drawBox([w_tot - w_leg, 0, w_leg], [w_leg, a, h_w], core_color);
    end

    % --- Center Post ---
    if strcmp(shape, 'round')
        [X, Y, Z] = cylinder(b/2, 50);
        Z = Z * leg_top_z;
        X = X + x_center;
        Y = Y + y_center;
        surf(X, Y, Z, 'FaceColor', core_color, 'EdgeColor', 'none');
    else
        drawBox([x_center - b/2, 0, w_leg], [b, a, h_w - lg], core_color);
    end

    % --- Winding ---
    z_bottom = w_leg;
    z_top    = z_bottom + h_pcb;

    if strcmp(shape, 'round')
        % Annular winding around round center post
        r_inner  = b/2 + s_ct;
        r_outer  = r_inner + w_winding;
        n_theta  = 60;
        theta    = linspace(0, 2*pi, n_theta + 1);
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
        % Outer wall
        [X_cyl, Y_cyl, Z_cyl] = cylinder(r_outer, n_theta);
        Z_cyl = Z_cyl * h_pcb + z_bottom;
        surf(X_cyl + x_center, Y_cyl + y_center, Z_cyl, ...
            'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);
        % Inner wall
        [X_cyl_in, Y_cyl_in, Z_cyl_in] = cylinder(r_inner, n_theta);
        Z_cyl_in = Z_cyl_in * h_pcb + z_bottom;
        surf(X_cyl_in + x_center, Y_cyl_in + y_center, Z_cyl_in, ...
            'FaceColor', copper_pri, 'EdgeColor', 'k', 'LineWidth', 0.3);

    else
        % Rectangular winding around rectangular center post
        % Inner rectangle: post face + clearance on all sides
        x_inn = x_center - b/2 - s_ct;
        y_inn = 0        - s_ct;            % post starts at y=0
        x_out = x_center - b/2 - s_ct - w_winding;
        y_out = 0        - s_ct - w_winding;

        % Inner dims
        dx_inn = b + 2*s_ct;
        dy_inn = a + 2*s_ct;
        % Outer dims
        dx_out = b + 2*s_ct + 2*w_winding;
        dy_out = a + 2*s_ct + 2*w_winding;

        % Origin of inner and outer rectangles (bottom-left corner)
        xy_inn = [x_center - dx_inn/2,  y_center - dy_inn/2];
        xy_out = [x_center - dx_out/2,  y_center - dy_out/2];

        drawRectTube(xy_inn, dx_inn, dy_inn, ...
                     xy_out, dx_out, dy_out, ...
                     z_bottom, z_top, copper_pri);
    end

    % --- Formatting ---
    xlabel('X [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Z [mm]', 'FontSize', 12, 'FontWeight', 'bold');

    title_str = sprintf(['                           3D EQ-Core Inductor \n' ...
        '      Footprint Area = %.1f in^2,  Total Loss = %.2f W\n'], ...
        result.A_footprint * 1e6 * 0.00155, ...
        result.P_total);
    title(title_str, 'FontSize', 13, 'FontWeight', 'bold');

    % --- Info textbox ---
    info_str = sprintf([ ...
        '-- LOSSES --\n' ...
        'P_copper  = %.3f W\n' ...
        'P_core    = %.3f W\n' ...
        'P_total   = %.3f W\n' ...
        'Temperature = %.3f C\n' ...
        '\n-- GEOMETRY --\n' ...
        'Number Turns = %.0f\n' ...
        'w_tot [X]    = %.2f mm\n' ...
        'w_winding    = %.2f mm\n' ...
        'l_tot [Y]    = %.2f mm\n' ...
        'l_core       = %.2f mm\n' ...
        'h_core [Z]   = %.2f mm\n' ...
        'l_gap        = %.3f mm\n' ...
        'gap loc      = %s'], ...
        result.P_cond, ...
        result.P_core, ...
        result.P_total, ...
        result.T_Lr, ...
        result.N, ...
        result.w_tot     * 1e3, ...
        result.w_winding * 1e3, ...
        result.l         * 1e3, ...
        result.a         * 1e3, ...
        result.h         * 1e3, ...
        result.lg        * 1e3, ...
        gap_loc);

    info_panel = uipanel('Units', 'normalized', ...
        'Position',        [0.77, 0.25, 0.2, 0.45], ...
        'BackgroundColor', 'w', ...
        'BorderType',      'line', ...
        'HighlightColor',  'k', ...
        'FontSize',        10);
    uicontrol(info_panel, ...
        'Style',           'text', ...
        'Units',           'normalized', ...
        'Position',        [0.03, 0.02, 0.96, 0.96], ...
        'String',          info_str, ...
        'FontSize',        10, ...
        'FontName',        'FixedWidth', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w');

    rotate3d on;
    hold off;
end

% ═════════════════════════════════════════════════════════════════════════
% HELPER: draw a hollow rectangular tube (winding cross-section)
%   xy_inn / xy_out : [x0, y0] bottom-left corner of inner / outer rect
%   dx, dy          : width and depth of each rectangle
%   z0, z1          : bottom and top z
% ═════════════════════════════════════════════════════════════════════════
function drawRectTube(xy_inn, dx_inn, dy_inn, xy_out, dx_out, dy_out, z0, z1, col)

    xi = xy_inn(1);  yi = xy_inn(2);
    xo = xy_out(1);  yo = xy_out(2);

    % Corner coordinates of inner and outer rectangles
    inn = [xi,           yi;           % bottom-left
           xi + dx_inn,  yi;           % bottom-right
           xi + dx_inn,  yi + dy_inn;  % top-right
           xi,           yi + dy_inn]; % top-left

    out = [xo,           yo;
           xo + dx_out,  yo;
           xo + dx_out,  yo + dy_out;
           xo,           yo + dy_out];

    face_col = col;
    edge_col = 'k';
    lw       = 0.3;

    % Top and bottom annular faces (4 trapezoid panels each)
    for k = 1:4
        k2 = mod(k, 4) + 1;
        % Vertices of one panel: inner edge, then outer edge
        Xf = [inn(k,1), inn(k2,1), out(k2,1), out(k,1)];
        Yf = [inn(k,2), inn(k2,2), out(k2,2), out(k,2)];
        % Bottom face
        fill3(Xf, Yf, ones(1,4)*z0, face_col, 'EdgeColor', edge_col, 'LineWidth', lw);
        % Top face
        fill3(Xf, Yf, ones(1,4)*z1, face_col, 'EdgeColor', edge_col, 'LineWidth', lw);
    end

    % Four outer vertical walls
    for k = 1:4
        k2 = mod(k, 4) + 1;
        Xf = [out(k,1), out(k2,1), out(k2,1), out(k,1)];
        Yf = [out(k,2), out(k2,2), out(k2,2), out(k,2)];
        Zf = [z0, z0, z1, z1];
        fill3(Xf, Yf, Zf, face_col, 'EdgeColor', edge_col, 'LineWidth', lw);
    end

    % Four inner vertical walls
    for k = 1:4
        k2 = mod(k, 4) + 1;
        Xf = [inn(k,1), inn(k2,1), inn(k2,1), inn(k,1)];
        Yf = [inn(k,2), inn(k2,2), inn(k2,2), inn(k,2)];
        Zf = [z0, z0, z1, z1];
        fill3(Xf, Yf, Zf, face_col, 'EdgeColor', edge_col, 'LineWidth', lw);
    end
end