function VIRT3Dfigure_raceTrack(result, material, T_tx, idx)

    figure(1)
    if idx == 0
    elseif idx == 1; subplot(2,2,1);
    elseif idx == 2; subplot(2,2,2);
    elseif idx == 3; subplot(2,2,3);
    elseif idx == 4; subplot(2,2,4);
    end
    cla; hold on; axis equal; view(35,15); grid on; grid minor; box on;

    % --- Colors ---
    core_color = [0.5, 0.5, 0.8];
    copper_pri = [0.85, 0.45, 0.01];

    % --- Unpack geometry (mm) ---
    a         = result.a         * 1e3;
    b         = result.b         * 1e3;
    w         = result.w         * 1e3;
    s_ct      = result.s_ct      * 1e3;
    h_pcb     = result.h_pcb     * 1e3;
    lg        = result.lg        * 1e3;
    w_winding = result.w_winding * 1e3;
    h_w       = result.h_w       * 1e3;
    w_tot     = result.w_tot     * 1e3;
    h         = result.h         * 1e3;
    w_leg     = result.w_leg     * 1e3;
    gap_loc   = result.gap_loc;

    % x_center: same as round version
    x_center  = w_tot / 2;
    % Center post is centered at (x_center, a/2) — same as round
    y_center  = a / 2;

    leg_top_z = w_leg + h_w - lg;   % bottom of air gap, same as round

    % =========================================================
    % HELPER: draw a racetrack prism (solid, extruded in z)
    %   cx, cy = center
    %   rlen   = b/2 = semicircle radius (short axis)
    %   rstra  = (a - b)/2 = half straight-section length (long axis offset)
    %   z0, z1 = bottom and top
    %   col    = face color
    %
    %   Shape: rectangle (a x b) with semicircles of radius b/2 on the
    %   short (y) ends. Total extents: x in [cx-b/2, cx+b/2],
    %                                  y in [cy-a/2, cy+a/2]
    % =========================================================
    function drawRaceTrackPrism(cx, cy, rlen, rstra, z0, z1, col)
        n_arc   = 40;
        % Right semicircle cap (y+ end), centered at (cx, cy+rstra)
        theta_r = linspace(0, pi, n_arc);
        % Left semicircle cap (y- end), centered at (cx, cy-rstra)
        theta_l = linspace(pi, 2*pi, n_arc);

        xc = [ rlen*cos(theta_r) + cx,  rlen*cos(theta_l) + cx  ];
        yc = [ rlen*sin(theta_r) + (cy + rstra), ...
               rlen*sin(theta_l) + (cy - rstra) ];
        nc = numel(xc);

        % Bottom face
        fill3(xc, yc, z0*ones(1,nc), col, 'EdgeColor','none','FaceAlpha',0.9);
        % Top face
        fill3(xc, yc, z1*ones(1,nc), col, 'EdgeColor','none','FaceAlpha',0.9);
        % Side walls
        for k = 1:nc
            k2 = mod(k,nc)+1;
            fill3([xc(k),xc(k2),xc(k2),xc(k)], ...
                  [yc(k),yc(k2),yc(k2),yc(k)], ...
                  [z0,z0,z1,z1], col, 'EdgeColor','none','FaceAlpha',0.9);
        end
    end

    % =========================================================
    % HELPER: draw racetrack winding annulus (hollow shell)
    %   Same coordinate convention as above, inner and outer contours
    % =========================================================
    function drawRaceTrackWinding(cx, cy, rlen_in, rstra_in, rlen_out, rstra_out, z0, z1, col)
        n_arc   = 40;
        theta_r = linspace(0,  pi, n_arc);
        theta_l = linspace(pi, 2*pi, n_arc);

        % Outer contour
        xo = [ rlen_out*cos(theta_r)+cx,  rlen_out*cos(theta_l)+cx  ];
        yo = [ rlen_out*sin(theta_r)+(cy+rstra_out), ...
               rlen_out*sin(theta_l)+(cy-rstra_out) ];
        % Inner contour
        xi = [ rlen_in*cos(theta_r)+cx,   rlen_in*cos(theta_l)+cx   ];
        yi = [ rlen_in*sin(theta_r)+(cy+rstra_in), ...
               rlen_in*sin(theta_l)+(cy-rstra_in)  ];
        nc = numel(xo);

        % Annular top and bottom faces
        xf = [xo, fliplr(xi)];
        yf = [yo, fliplr(yi)];
        fill3(xf, yf, z1*ones(1,2*nc), col, 'EdgeColor','k','LineWidth',0.3,'FaceAlpha',0.95);
        fill3(xf, yf, z0*ones(1,2*nc), col, 'EdgeColor','k','LineWidth',0.3,'FaceAlpha',0.95);

        % Outer wall
        for k = 1:nc
            k2 = mod(k,nc)+1;
            fill3([xo(k),xo(k2),xo(k2),xo(k)], ...
                  [yo(k),yo(k2),yo(k2),yo(k)], ...
                  [z0,z0,z1,z1], col, 'EdgeColor','k','LineWidth',0.3,'FaceAlpha',0.95);
        end
        % Inner wall
        for k = 1:nc
            k2 = mod(k,nc)+1;
            fill3([xi(k),xi(k2),xi(k2),xi(k)], ...
                  [yi(k),yi(k2),yi(k2),yi(k)], ...
                  [z0,z0,z1,z1], col, 'EdgeColor','k','LineWidth',0.3,'FaceAlpha',0.95);
        end
    end

    % =========================================================
    % CORE — identical structure to round version
    % =========================================================

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

    % --- Center Post (racetrack, gapped) ---
    % rlen  = b/2  (semicircle radius, x-direction)
    % rstra = (a - b)/2  (half straight offset, y-direction)
    rlen  = b / 2;
    rstra = (a - b) / 2;
    drawRaceTrackPrism(x_center, y_center, rlen, rstra, w_leg, leg_top_z, core_color);

    % =========================================================
    % WINDING — racetrack annulus around center post
    % The straight section offset (rstra) stays fixed — the winding
    % only grows outward in the radial (x) direction.
    % rlen grows by s_ct then w_winding; rstra stays at (a-b)/2.
    % =========================================================
    rlen_in   = rlen  + s_ct;
    rstra_in  = rstra;              % same straight offset as post
    rlen_out  = rlen_in + w_winding;
    rstra_out = rstra;              % same straight offset as post

    z_wind_bot = w_leg;
    z_wind_top = w_leg + h_pcb;

    drawRaceTrackWinding(x_center, y_center, ...
        rlen_in,  rstra_in,  ...
        rlen_out, rstra_out, ...
        z_wind_bot, z_wind_top, copper_pri);

    % --- Labels ---
    xlabel('X [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Z [mm]', 'FontSize', 12, 'FontWeight', 'bold');

    title_str = sprintf(['                     3D Racetrack VIRT Transformer\n' ...
        '      Footprint = %.1f in²,  Total Loss = %.2f W\n' ...
        '                              Temp = %.2f °C'], ...
        result.A_footprint * 1e6 * 0.00155, ...
        result.P_total, T_tx);
    title(title_str, 'FontSize', 13, 'FontWeight', 'bold');

    % --- Info textbox ---
    info_str = sprintf([ ...
        '-- LOSSES --\n' ...
        'P_pri    = %.3f W\n' ...
        'P_sec    = %.3f W\n' ...
        'P_core   = %.3f W\n' ...
        'P_total  = %.3f W\n' ...
        '\n-- GEOMETRY --\n' ...
        'a [inner len]  = %.2f mm\n' ...
        'b [inner dia]  = %.2f mm\n' ...
        'w [wall]       = %.2f mm\n' ...
        'w_tot [X]      = %.2f mm\n' ...
        'l_tot [Y]      = %.2f mm\n' ...
        'h_core [Z]     = %.2f mm\n' ...
        'l_gap          = %.3f mm\n' ...
        'gap loc        = %s'], ...
        result.P_pri,      result.P_sec, ...
        result.P_core,     result.P_total, ...
        result.a  * 1e3,   result.b  * 1e3, ...
        result.w  * 1e3,   result.w_tot * 1e3, ...
        result.l  * 1e3,   result.h  * 1e3, ...
        result.lg * 1e3,   gap_loc);

    info_panel = uipanel('Units','normalized', ...
        'Position',       [0.77, 0.25, 0.2, 0.45], ...
        'BackgroundColor','w', 'BorderType','line', ...
        'HighlightColor', 'k', 'FontSize', 10);
    uicontrol(info_panel, 'Style','text', 'Units','normalized', ...
        'Position',           [0.03, 0.02, 0.96, 0.96], ...
        'String',             info_str, 'FontSize', 10, ...
        'FontName',           'FixedWidth', ...
        'HorizontalAlignment','left', 'BackgroundColor','w');

    rotate3d on;
    hold off;
end