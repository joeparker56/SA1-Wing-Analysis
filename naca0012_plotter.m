%    foil.m
%
%  Script to analyse an aerofoil section using potential flow calculation
%  and boundary layer solver.
%

clear
close all

global Re

%  Read in the parameter file
caseref = 'naca0012';
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
alpha = linspace(-14    ,24,40);
cls = zeros(1,12);
cds = zeros(1,12);
lds = zeros(1,12);
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
   cls(nalpha) = Cl;
   cds(nalpha) = Cd;
   lds(nalpha) = Cl/Cd;

%    copy Cl and Cd into arrays for alpha sweep plots

  clswp(nalpha) = Cl;
  cdswp(nalpha) = Cd;
  lovdswp(nalpha) = Cl/Cd;

end

% 1. Load and display the uncropped image
% (imshow is preferred here over imagesc to maintain the exact aspect ratio)
img = imread('NACA0012_CL_ALPHA.png');
figure;
imshow(img);
hold on;

% 2. Calibrate the X and Y axes interactively
% The command window will prompt you to click the image.
disp('Click exactly on the BOTTOM-LEFT intersection of the graph axes.');
[px1, py1] = ginput(1);

disp('Click exactly on the TOP-RIGHT corner/end of the graph axes.');
[px2, py2] = ginput(1);

% Optional: Plot temporary markers so you can see where you clicked

% 3. Define the real-world data values for the points you just clicked
dx1 = -24;   % X-axis value at the bottom-left click
dy1 = -1.8;   % Y-axis value at the bottom-left click

dx2 = 32; % X-axis value at the top-right click
dy2 = 2.8;  % Y-axis value at the top-right click

% 4. Transform your MATLAB data into pixel coordinates
% (Assuming you already have 'x_data' and 'y_data' in your workspace)
% This formula automatically handles the fact that image pixels count from the top-down.
x_data = alpha;
y_data = cls;
x_mapped = px1 + (x_data - dx1) * (px2 - px1) / (dx2 - dx1);
y_mapped = py1 + (y_data - dy1) * (py2 - py1) / (dy2 - dy1);

% 5. Plot your mapped data directly over the image
plot(x_mapped, y_mapped, 'rx', 'MarkerSize', 4, 'LineWidth',1);

% Optional: Add a legend
legend('Panel method solver for Re=3e6');
hold off;
print -deps2c Joe_Report/NACA0012_CL_A.eps;




img = imread('NACA0012_CD_CL.png');
figure(2);
imshow(img);
hold on;

% 2. Calibrate the X and Y axes interactively
% The command window will prompt you to click the image.
disp('Click exactly on the BOTTOM-LEFT intersection of the graph axes.');
[px1, py1] = ginput(1);

disp('Click exactly on the TOP-RIGHT corner/end of the graph axes.');
[px2, py2] = ginput(1);

% Optional: Plot temporary markers so you can see where you clicked

% 3. Define the real-world data values for the points you just clicked
dx1 = -1.6;   % X-axis value at the bottom-left click
dy1 = -0.014;   % Y-axis value at the bottom-left click

dx2 = 1.6; % X-axis value at the top-right click
dy2 = 0.032;  % Y-axis value at the top-right click

% 4. Transform your MATLAB data into pixel coordinates
% (Assuming you already have 'x_data' and 'y_data' in your workspace)
% This formula automatically handles the fact that image pixels count from the top-down.
x_data = cls;
y_data = cds;
x_mapped = px1 + (x_data - dx1) * (px2 - px1) / (dx2 - dx1);
y_mapped = py1 + (y_data - dy1) * (py2 - py1) / (dy2 - dy1);

% 5. Plot your mapped data directly over the image
plot(x_mapped, y_mapped, 'rx', 'MarkerSize', 4, 'LineWidth',1);

% Optional: Add a legend
legend('Panel method solver for Re=3e6');
hold off;
print -deps2c Joe_Report/NACA0012_CD_CL.eps;
