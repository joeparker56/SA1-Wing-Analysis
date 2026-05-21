clear;
close all;


global re_l dudx u0

re_l = 1e7;
dudx =  -0;
u0 = 1;

function dthickdx = thickdash(xmx0, thick)
    global re_l dudx u0
    dthickdx = zeros(2,1);
    He = thick(2)/thick(1);
    if He>=1.46
        H =(11*He + 15)/(48*He-59);
    else
        H = 2.803;
    end 
    ue = u0 + dudx*xmx0;
    re_thet = re_l * ue * thick(1);
    cf = 0.091416*((H-1)*re_thet)^(-0.232) * exp(-1.26*H);
    cdiss = 0.010026*((H-1)*re_thet)^(-1/6);
    dthickdx(1,1) = cf/2 - ((H+2)/ue)*dudx*thick(1);
    dthickdx(2,1) = cdiss - (3/ue)*dudx*thick(2);
end

x0 = 0.01;
thick0(1) = 0.023*x0*(re_l*x0)^(-1/6);
thick0(2) = 1.83*thick0(1);

[delx, thickhist] =ode45(@thickdash, [0, 0.99], thick0);

theta7 = 0.037*delx.*(re_l.*(delx+x0)).^(-1/5);
theta9 = 0.023*delx.*(re_l.*(delx+x0)).^(-1/6);

figure;
plot(delx+x0, thickhist(:,1))
hold on;
plot(delx+x0, theta7)
hold off;
hold on;
plot(delx+x0, theta9)
hold off;
