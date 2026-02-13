function drawBox(origin, dimensions, color)
% Draw a solid rectangular box as a patch object.
%
% Inputs:
%   origin     - [x0, y0, z0] position of the box corner [mm]
%   dimensions - [dx, dy, dz] side lengths of the box [mm]
%   color      - [r, g, b] face color

    x0 = origin(1);     y0 = origin(2);     z0 = origin(3);
    dx = dimensions(1); dy = dimensions(2); dz = dimensions(3);

    vertices = [
        x0,    y0,    z0;
        x0+dx, y0,    z0;
        x0+dx, y0+dy, z0;
        x0,    y0+dy, z0;
        x0,    y0,    z0+dz;
        x0+dx, y0,    z0+dz;
        x0+dx, y0+dy, z0+dz;
        x0,    y0+dy, z0+dz;
    ];

    faces = [
        1 2 6 5;    % Front
        2 3 7 6;    % Right
        3 4 8 7;    % Back
        4 1 5 8;    % Left
        1 2 3 4;    % Bottom
        5 6 7 8;    % Top
    ];

    patch('Vertices', vertices, 'Faces', faces, ...
          'FaceColor', color, 'EdgeColor', 'k', ...
          'LineWidth', 1, 'FaceAlpha', 0.9);

end