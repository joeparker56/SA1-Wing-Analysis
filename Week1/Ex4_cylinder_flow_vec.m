clear;
close all;

xmin = -2.5; xmax = 2.5;
ymin = -2;   ymax = 2;
ny = 41; nx = 51;

np = 100;
theta = (0:np)*2*pi/np;
xs = cos(theta);
ys = sin(theta);
gamma_panel_ends = -2*sin(theta);

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
[xm, ym] = meshgrid(x, y);

psi = ym;

for s = 1:np
    [infa, infb] = panelinf_vec(xs(s), ys(s), xs(s+1), ys(s+1), xm, ym);
    psi = psi + gamma_panel_ends(s)*infa + gamma_panel_ends(s+1)*infb;
end

c = -1.75:0.25:1.75;
contour(xm, ym, psi, c);
xlabel('x'); ylabel('y');
title('Contour plot of flow around cylinder');
hold on
plot(xs, ys)
hold off