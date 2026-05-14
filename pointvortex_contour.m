clear;
close all;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
ny = 41;
nx = 51;
Gamma = 3.0;
xc = 0.5;
yc = 0.25;

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(xmax-xmin)/(ny-1);
        psi(i,j) = psipv(xc, yc, Gamma, xm(i,j), ym(i,j));
    end
end

c = -0.4:0.2:1.2;
contour(xm,ym,psi,c);
