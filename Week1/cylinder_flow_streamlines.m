

xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
ny = 41;
nx = 51;
del = 1.5;

psi = zeros(nx, ny);
xm = zeros(nx, ny);
ym = zeros(nx, ny);


np = 100;
theta = (0:np)*2*pi/np;

for i = 1:np+1
    xs(i) = cos(theta(i));
    ys(i) = sin(theta(i));
end


gamma_panel_ends = -2 *sin(theta);

disp(gamma_panel_ends)

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        psi(i,j) = ym(i,j);
        for s = 1:np
            [infa, infb] = panelinf(xs(s), ys(s), xs(s+1), ys(s+1), xm(i,j), ym(i,j));
            psi(i,j) = psi(i,j) + gamma_panel_ends(s)*infa + gamma_panel_ends(s+1) *infb;
        end
    end
end


c = -1.75:0.25:1.75;


contour(xm,ym,psi,c);
xlabel('x');
ylabel('y');
title('Contour plot of flow around cylinder');

hold on
plot(xs,ys)
hold off