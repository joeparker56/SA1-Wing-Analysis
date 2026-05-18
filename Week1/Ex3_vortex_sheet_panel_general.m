clear;
close all;
xmin = 0;
xmax = 5;
ymin = 0;
ymax = 4;
ny = 41;
nx = 51;
del = 1.5;
%gamma_a = 0.5;
%gamma_b = 1.5;
nv=100;
yc=0;
psi_approx = zeros(nx, ny);
infa = zeros(nx,ny) ;
infb = zeros(nx, ny);
xm = zeros(nx, ny);
ym = zeros(nx, ny);
psi_analy = zeros(nx,ny);
infa_approx = zeros(nx,ny);
infb_approx = zeros(nx,ny);
xa = 3.5;
ya = 2.5;
xb = 1.6;
yb = 1.1;

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        [infa(i,j), infb(i,j)] = panelinf(xa, ya, xb, yb, xm(i,j), ym(i,j));
        %psi_analy(i,j) = gamma_a*infa(i,j) + gamma_b *infb(i,j); 
      
    end
end


l=0;
panel_length = norm([(xb-xa), (yb-ya)]);
dl = panel_length/nv; 
pl_vec = [(xb-xa), (yb-ya)]/norm([(xb-xa), (yb-ya)]);
a = [xa,ya];
b = [xb, yb];
gamma_a = 1.0;
gamma_b = 0;

%evaluate approximate influence coefficient a
for n = 1:nv
    %this code is designed to evaluate arbitrary values of gamma_a and
    %gamma_b but has been adapted to just evaluate influence a and
    %influence b
    centre = a + ((2*n-1)/2)*dl*pl_vec;
    Gamma = (gamma_a*norm(b-centre)/panel_length + gamma_b*norm(centre-a)/panel_length)*dl;
    for i = 1:nx
        for j = 1:ny
            infa_approx(i,j) = infa_approx(i,j) + psipv(centre(1), centre(2), Gamma, xm(i,j), ym(i,j));
        end 
    end
end

%evaluate approximate influence coefficient b
gamma_b = 1.0;
gamma_a = 0;
for n = 1:nv
    centre = a + ((2*n-1)/2)*dl*pl_vec;
    Gamma = (gamma_a*norm(b-centre)/panel_length + gamma_b*norm(centre-a)/panel_length)*dl;
    for i = 1:nx
        for j = 1:ny
            infb_approx(i,j) = infb_approx(i,j) + psipv(centre(1), centre(2), Gamma, xm(i,j), ym(i,j));
        end 
    end
end

c = -0.15:0.05:0.15;





contour(xm,ym,infa,c);
xlabel('x');
ylabel('y');
title('Contour plot of influence coefficient a');
print -deps2c infa_ex3.eps;


figure(2);
contour(xm,ym,infb,c);
xlabel('x');
ylabel('y');
title('Contour plot of influence coefficient b');
print -deps2c infb_ex3.eps;


figure(3);
contour(xm,ym,infa_approx,c);
xlabel('x');
ylabel('y');
title('Contour plot of the approximate influence coefficient a');
print -deps2c infa_approx_ex3.eps;


figure(4);
contour(xm,ym,infb_approx,c);
xlabel('x');
ylabel('y');
title('Contour plot of the approximate influence coefficient b');
print -deps2c infb_approx_ex3.eps;
