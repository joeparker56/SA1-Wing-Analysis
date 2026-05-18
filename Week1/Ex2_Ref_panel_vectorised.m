xmin = -2.5; xmax = 2.5;
ymin = -2;   ymax = 2;
ny = 41; nx = 51;
del = 1.5;
gamma_a = 0.5; gamma_b = 1.5;
nv = 100;
yc = 0;

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
[xm, ym] = meshgrid(x, y);

[infa, infb] = refpainf_vec(del, xm, ym);
psi_analy = gamma_a*infa + gamma_b*infb;

psi_approx = zeros(size(xm));
   
for n = 1:nv
    xc = (2*n-1)*del/(2*nv);
    Gamma = (gamma_a*(del-xc)/del + gamma_b*xc/del)*del/nv;
    psi_approx = psi_approx + psipv(xc, yc, Gamma, xm, ym);
end



contour(xm,ym,infa,c);
xlabel('x');
ylabel('y');
title('Contour plot of influence coefficient a');

figure(2);
contour(xm,ym,infb,c);
xlabel('x');
ylabel('y');
title('Contour plot of influence coefficient b');

figure(3);
contour(xm,ym,psi_analy,c);
xlabel('x');
ylabel('y');
title('Contour plot of the analytical solution');

figure(4);
contour(xm,ym,psi_approx,c);
xlabel('x');
ylabel('y');
title('Contour plot of discretised approximation');