clear;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
ny = 41;
nx = 51;

psi = zeros(nx, ny);
psi_fs = zeros(nx,ny);
xm = zeros(nx, ny);
ym = zeros(nx, ny);
alpha = pi/24;

np = 100;
theta = (0:np)*2*pi/np;
xs = zeros(np+1);
ys = zeros(np+1);

for i = 1:np+1
    xs(i) = cos(theta(i));
    ys(i) = sin(theta(i));
end


A = build_lhs(xs, ys);
b = build_rhs(xs, ys, alpha);
gam = A\b;


for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        psi(i,j) = ym(i,j)*cos(alpha) - xm(i,j)*sin(alpha); %free stream steamfunction
        for s = 1:np
            [infa, infb] = panelinf(xs(s), ys(s), xs(s+1), ys(s+1), xm(i,j), ym(i,j));
            psi(i,j) = psi(i,j) + gam(s)*infa + gam(s+1) *infb;
        end
    end
end


c = -1.75:0.25:1.75;
circ = 0;
for i=1:np
    circ = circ + gam(i)*2*pi/np;
end

contour(xm,ym,psi,c);
xlabel('x');
ylabel('y');
title('Contour plot of flow around cylinder for alpha = pi/24');

hold on
plot(xs,ys)
hold off
print -deps2c contour_alpha_pi24.eps;


theta_pi = theta/pi;
figure(2);
plot(theta_pi,gam);
xlabel('θ/π');
ylabel('Surface velocity (γ)');
title('Plot of surface velocity against theta for alpha = pi/24')
axis([0 2 -2.5 2.5]);
print -deps2c sur_vel_alpha_pi24.eps;
