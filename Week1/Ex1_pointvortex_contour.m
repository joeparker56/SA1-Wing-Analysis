clear;
close all;

xmin = -2.5; xmax = 2.5;
ymin = -2;   ymax = 2;
nx = 51; ny = 41;
Gamma = 3.0;
xc = 0.5; yc = 0.25;

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);

[xm, ym] = meshgrid(x, y);

psi = psipv(xc, yc, Gamma, xm, ym);

c = -0.4:0.2:1.2;
contour(xm, ym, psi, c);
print -deps2c point_vortex_contour_ex1.eps;
