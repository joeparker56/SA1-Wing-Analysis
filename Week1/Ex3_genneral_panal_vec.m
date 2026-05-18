clear;
close all;

xmin = 0; xmax = 5;
ymin = 0; ymax = 4;
ny = 41; nx = 51;
del = 1.5;
gamma_a = 0.5; gamma_b = 1.5;
nv = 100;
yc = 0;

xa = 3.5; ya = 2.5;
xb = 1.6; yb = 1.1;

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
[xm, ym] = meshgrid(x, y);

[infa, infb] = panelinf_vec(xa, ya, xb, yb, xm, ym);
psi_analy = gamma_a*infa + gamma_b*infb;

psi_approx = zeros(size(xm));

panel_length = norm([xb-xa, yb-ya]);
dl = panel_length/nv;
pl_vec = [xb-xa, yb-ya]/panel_length;
a = [xa, ya];
b = [xb, yb];

for n = 1:nv
    centre = a + ((2*n-1)/2)*dl*pl_vec;
    Gamma = (gamma_a*norm(b-centre)/panel_length + gamma_b*norm(centre-a)/panel_length)*dl;
    psi_approx = psi_approx + psipv(centre(1), centre(2), Gamma, xm, ym);
end

c = -0.15:0.05:0.15;

contour(xm, ym, infa, c);
xlabel('x'); ylabel('y');
title('Contour plot of influence coefficient a');

figure(2);
contour(xm, ym, infb, c);
xlabel('x'); ylabel('y');
title('Contour plot of influence coefficient b');

figure(3);
contour(xm, ym, psi_analy, c);
xlabel('x'); ylabel('y');
title('Contour plot of the analytical solution');

figure(4);
contour(xm, ym, psi_approx, c);
xlabel('x'); ylabel('y');
title('Contour plot of discretised approximation');