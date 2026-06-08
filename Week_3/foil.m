%    foil.m
%
%  Script to analyse an aerofoil section using potential flow calculation
%  and boundary layer solver.
%

clear
close all

global Re

%  Read in the parameter file
caseref = input('\nEnter case reference: ','s');
parfile = ['Parfiles/' caseref '.txt'];
fprintf(1, '%s\n\n', ['Reading in parameter file: ' parfile])
[section np Re alpha] = par_read(parfile);

%  Read in the section geometry
secfile = ['Geometry/' section '.surf'];
[xk yk] = textread ( secfile, '%f%f' );

%  Generate high-resolution surface description via cubic splines
nphr = 5*np;
[xshr yshr] = splinefit ( xk, yk, nphr );

%  Resize section so that it lies between (0,0) and (1,0)
[xsin ysin] = resyze ( xshr, yshr );

%  Interpolate to required number of panels (uniform size)
[xs ys] = make_upanels ( xsin, ysin, np );

%  Assemble the lhs of the equations for the potential flow calculation
A = build_lhs_fully_vec ( xs, ys );
Am1 = inv(A);

%  Loop over alpha values
for nalpha = 1:length(alpha)
    
%    rhs of equations
  alfrad = pi * alpha(nalpha)/180;
  b = build_rhs_vec ( xs, ys, alfrad );

%    solve for surface vortex sheet strength
  gam = Am1 * b;

%    calculate cp distribution and overall circulation
  [cp circ] = potential_op ( xs, ys, gam );

  
%    locate stagnation point and calculate stagnation panel length
  [ipstag fracstag] = find_stag(gam);
  dsstag = sqrt((xs(ipstag+1)-xs(ipstag))^2 + (ys(ipstag+1)-ys(ipstag))^2);

%    upper surface boundary layer calc

%    first assemble pressure distribution along bl
  clear su cpu
  su(1) = fracstag*dsstag;
  cpu(1) = cp(ipstag);
  for is = ipstag-1:-1:1
    iu = ipstag - is + 1;
    su(iu) = su(iu-1) + sqrt((xs(is+1)-xs(is))^2 + (ys(is+1)-ys(is))^2);
    cpu(iu) = cp(is);
  end

%    check for stagnation point at end of stagnation panel
  if fracstag < 1e-6
    su(1) = 0.01*su(2);    % go just downstream of stagnation
    uejds = 0.01 * sqrt(1-cpu(2));
    cpu(1) = 1 - uejds^2;
  end

%    boundary layer solver
  [iunt iuls iutr iuts delstaru thetau] = bl_solv ( su, cpu );

  %    lower surface boundary layer calc

%    first assemble pressure distribution along bl
  clear sl cpl
  sl(1) = (1-fracstag) * dsstag;
  cpl(1) = cp(ipstag+1);
  for is = ipstag+2:np+1
    il = is - ipstag;
    sl(il) = sl(il-1) + sqrt((xs(is-1)-xs(is))^2 + (ys(is-1)-ys(is))^2);
    cpl(il) = cp(is);
  end

%    check for stagnation point at end of stagnation panel
  if fracstag > 0.999999
    sl(1) = 0.01*sl(2);    % go just downstream of stagnation
    uejds = 0.01 * sqrt(1-cpl(2));
    cpl(1) = 1 - uejds^2;
  end

%    boundary layer solver
  [ilnt ills iltr ilts delstarl thetal] = bl_solv ( sl, cpl );
    

%    lift and drag coefficients
  [Cl Cd] = forces ( circ, cp, delstarl, thetal, delstaru, thetau );

%    copy Cl and Cd into arrays for alpha sweep plots

  clswp(nalpha) = Cl;
  cdswp(nalpha) = Cd;
  lovdswp(nalpha) = Cl/Cd;

%    screen output

  disp ( sprintf ( '\n%s%5.3f%s', ...
                   'Results for alpha = ', alpha(nalpha), ' degrees' ) )

  disp ( sprintf ( '\n%s%5.3f', '  Lift coefficient: ', Cl ) )
  disp ( sprintf ( '%s%7.5f', '  Drag coefficient: ', Cd ) )
  disp ( sprintf ( '%s%5.3f\n', '  Lift-to-drag ratio: ', Cl/Cd ) )

  upperbl = sprintf ( '%s', '  Upper surface boundary layer:' );
  if iunt~=0
    is = ipstag + 1 - iunt;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Natural transition at x = ', xs(is) );
      iunt = ipstag + 1 - iunt;

  end
  if iuls~=0
    is = ipstag + 1 - iuls;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Laminar separation at x = ', xs(is) );
      iuls = ipstag + 1 - iuls;

    if iutr~=0
      is = ipstag + 1 - iutr;
      upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                          '    Turbulent reattachment at x = ', xs(is) );
        iutr = ipstag + 1 - iutr;

    end
  end
  if iuts~=0
    is = ipstag + 1 - iuts;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Turbulent separation at x = ', xs(is) );
      iuts = ipstag + 1 - iuts;

  end
        

  upperbl = sprintf ( '%s\n', upperbl );
  disp(upperbl)



  lowerbl = sprintf ( '%s', '  Lower surface boundary layer:' );
  if ilnt~=0
    is = ipstag + ilnt;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Natural transition at x = ', xs(is) );
      ilnt = ipstag + ilnt;

  end
  if ills~=0
    is = ipstag + ills;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Laminar separation at x = ', xs(is) );
      ills = ipstag + ills;
    if iltr~=0
      is = ipstag + iltr;
      lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                          '    Turbulent reattachment at x = ', xs(is) );
        iltr = ipstag + iltr;

    end
  end
  if ilts~=0
    is = ipstag + ilts;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Turbulent separation at x = ', xs(is) );
      ilts = ipstag + ilts;

  end
  fig = figure;
  
  % Create a stacked layout (3 rows total: top plot gets 2, bottom gets 1)
  tiledlayout(3, 1, 'TileSpacing', 'compact');
  
  % Define the global X-axis bounds so both plots match perfectly
  xmin = -0.4;
  xmax = 1.4;

  % ==========================================
  % TOP PLOT: -Cp Distribution
  % ==========================================
  ax1 = nexttile([2 1]);
  hold(ax1, 'on');
  
  % Suction side (upper surface) in blue, pressure side (lower) in red
  % (split at the stagnation panel; include ipstag in both so they join)
  plot(ax1, xs(1:ipstag),   -cp(1:ipstag),   'b-', 'LineWidth', 1.5)
  plot(ax1, xs(ipstag:end), -cp(ipstag:end), 'r-', 'LineWidth', 1.5)
  
  % Add boundary layer labels
  add_bl_labels(ax1, xs, -cp, iunt, iuls, iutr, iuts, 'blue', 'bottom');
  add_bl_labels(ax1, xs, -cp, ilnt, ills, iltr, ilts, 'red', 'top');
  
  % Data Box - left-aligned text, pinned to the top-right corner
  info_box = sprintf('C_L = %.4f\nC_D = %.6f\nL/D = %.2f', Cl, Cd, Cl/Cd);
  t = text(ax1, 0.95, 0.95, info_box, 'Units', 'normalized', ...
      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
      'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 5, 'FontSize', 10);
  % Shift left by the box width so its right edge stays in the corner
  drawnow;
  t.Position(1) = 0.95 - t.Extent(3);
      
  ylabel(ax1, '-Cp')
  title(ax1, ['alpha = ' num2str(alpha(nalpha))])
  grid(ax1, 'on');
  hold(ax1, 'off');

  % ==========================================
  % BOTTOM PLOT: Airfoil & Streamlines
  % ==========================================
  ax2 = nexttile;
  hold(ax2, 'on');
  
  % 1. Streamlines
  ymin = -max(ys)*1.5;
  ymax = -ymin;
  x_grid = linspace(xmin, xmax, 100);
  y_grid = linspace(ymin, ymax, 100);
  [xm, ym] = meshgrid(x_grid, y_grid);
  
  psi = ym*cos(alfrad) - xm*sin(alfrad);
  for s = 1:np
      [infa, infb] = panelinf_vec(xs(s), ys(s), xs(s+1), ys(s+1), xm, ym);
      psi = psi + gam(s)*infa + gam(s+1)*infb;
  end

  % Velocity field from the streamfunction: u = d(psi)/dy, v = -d(psi)/dx
  [dpsidx, dpsidy] = gradient(psi, x_grid, y_grid);
  u = dpsidy;
  v = -dpsidx;
  speed = sqrt(u.^2 + v.^2);

  % Blank out points inside the airfoil (psi there is non-physical)
  speed(inpolygon(xm, ym, xs, ys)) = NaN;

  % Colour-coded velocity-magnitude background.
  % Use fixed levels and colour limits so a given velocity (e.g. the
  % free-stream value of 1) always maps to the same colour across cases.
  vmax = 2;                              % top of the colour scale (V/Vinf)
  vlevels = linspace(0, vmax, 41);
  contourf(ax2, xm, ym, speed, vlevels, 'LineColor', 'none');
  colormap(ax2, 'jet');
  clim(ax2, [0 vmax]);                   % lock the colour axis (use caxis on older MATLAB)
  cb = colorbar(ax2);
  ylabel(cb, 'u/U_\infty');

  % Streamlines overlaid as thin black contours of psi
  c = -0.5:0.025:0.5;
  contour(ax2, xm, ym, psi, c, 'LineColor', 'k', 'LineWidth', 0.5);

  % 2. Airfoil Geometry
  plot(ax2, xs, ys, 'Color', 'black', 'LineWidth', 2)
    
  add_bl_labels(ax2, xs, ys, iunt, iuls, iutr, iuts, 'blue', 'bottom');
  add_bl_labels(ax2, xs, ys, ilnt, ills, iltr, ilts, 'red', 'top');
  % 3. Stagnation Point
  x_stag = xs(ipstag) + fracstag * (xs(ipstag+1) - xs(ipstag));
  y_stag = ys(ipstag) + fracstag * (ys(ipstag+1) - ys(ipstag));
  plot(ax2, x_stag, y_stag, 'o', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], 'MarkerSize', 6);
  text(ax2, x_stag - 0.03, y_stag, 'Stag', 'Color', [0 0.6 0], ...
      'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontSize', 9, ...
      'FontWeight', 'bold', 'BackgroundColor', 'white', 'EdgeColor', [0 0.6 0], 'Margin', 1);

  ylabel(ax2, 'y')
  xlabel(ax2, 'x')
  
  % Set physical scaling to 1:1 so the airfoil isn't squashed
  % Do NOT use "axis equal" here!
  set(ax2, 'DataAspectRatio', [1 1 1]); 
  ylim(ax2, [ymin ymax]);
  grid(ax2, 'on');
  hold(ax2, 'off');

  % ==========================================
  % ALIGNMENT & SAVING
  % ==========================================
  % Link the axes and force the exact same X limits on both
  linkaxes([ax1, ax2], 'x');
  xlim(ax1, [xmin xmax]); 

  % save data for this alpha
  fname = ['Data/' caseref '_' num2str(alpha(nalpha)) '.mat'];
  save(fname, 'Cl', 'Cd', 'xs', 'cp', ...
         'sl', 'delstarl', 'thetal', 'lowerbl', ...
         'su', 'delstaru', 'thetau', 'upperbl');
  location = ['Joe_Report/plots/' caseref '_' num2str(alpha(nalpha))];
  print(fig, location, '-deps2c');
end

%  save alpha sweep data in summary file

fname = ['Data/' caseref '.mat'];
save ( fname, 'xs', 'ys', 'alpha', 'clswp', 'cdswp', 'lovdswp' )

function add_bl_labels(ax,x, y, i_nt, i_ls, i_tr, i_ts, color, alignment)

    % Gather the active labels for this surface
    idx = [i_nt, i_ls, i_tr, i_ts];
    txt = {'NT', 'LS', 'TR', 'TS'};
    keep = idx > 0 & ~isnan(idx);
    idx = idx(keep);
    lab = txt(keep);
    if isempty(idx)
        return
    end

    % Modest, fixed vertical offset to lift labels just clear of the curve,
    % scaled to this plot's y-range so it works on both axes
    off = 0.05 * (max(y) - min(y));
    if ~strcmp(alignment, 'bottom')
        off = -off;   % below the marker for lower-surface (top alignment)
    end

    % Sort by x so neighbours that sit on top of each other can be detected
    xp = x(idx);  xp = xp(:);
    yp = y(idx);  yp = yp(:);
    [xp, order] = sort(xp);
    yp  = yp(order);
    lab = lab(order);

    % Group labels whose x positions nearly coincide, then fan each group out
    % sideways so the chips sit next to each other (short diagonal leaders)
    xthresh = 0.04 * (max(x) - min(x));
    dxstep  = 0.05 * (max(x) - min(x));
    grp = ones(numel(xp), 1);
    for k = 2:numel(xp)
        if abs(xp(k) - xp(k-1)) < xthresh
            grp(k) = grp(k-1);          % same cluster
        else
            grp(k) = grp(k-1) + 1;      % new cluster
        end
    end
    for g = unique(grp)'
        m = find(grp == g);
        n = numel(m);
        for j = 1:n
            hdx = (j - (n+1)/2) * dxstep;   % fan symmetrically about the marker
            place_label(ax, xp(m(j)), yp(m(j)), hdx, off, lab{m(j)}, color, alignment);
        end
    end
end

function place_label(ax, xp, yp, dx, dy, txt, color, alignment)
    % Marker on the surface
    plot(ax, xp, yp, 'o', 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 5);
    % Thin leader line from the marker to the offset label
    plot(ax, [xp xp+dx], [yp yp+dy], 'Color', color, 'LineWidth', 0.5);
    % Label in a white chip for constant contrast against the colour field
    text(ax, xp+dx, yp+dy, txt, 'Color', color, 'FontWeight', 'bold', ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', alignment, ...
        'BackgroundColor', 'white', 'EdgeColor', color, 'Margin', 1);
end