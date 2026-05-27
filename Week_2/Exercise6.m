clear;
close all;

global re_l dudx u0

np = 101;
%for turbulent seperation at x/L=1, dudx = -0.383
x = linspace(0,1,np); %length varies from 0->1
dudx = 0.; %dimensionless gradient



re_l = 1e7;

theta = zeros(1,np);  %theta squared at each position
Hes = zeros(1,np);
Hes(1) = 1.57258;


laminar = true;
integ = 0;
int = 0; %location of laminar->turbulent transition
ils = 0; %location of laminar seperation
itr = 0; %location of turbulent reattachment 
its = 0; %location of turbulent seperation
i=1;
u_ = linspace(1, 1+dudx, np);
while laminar && i<np
    i=i+1;
    integ = integ + ueintbit(x(i-1), u_(i-1), x(i), u_(i));
    theta(i) = sqrt((0.45/re_l) * u_(i)^(-6) * integ);
    m = -re_l * theta(i)^2  * dudx;
    H = thwaites_lookup(m);
    Hes(i) = laminar_He(H);
    re_thet = re_l * u_(i) * theta(i);
    if log(re_thet) >= 18.4*Hes(i) - 21.74
        laminar = false;    %transition
        int = i;

    elseif m>= 0.09
        laminar = false;
        ils = i;
        Hes(i) = 1.51509;
    end            
end
delta_e = Hes(i)*theta(i);

while its == 0 && i<np %loop for turbulent boundary layer
    i=i+1;
    thick0 = [theta(i-1), delta_e];
    u0 = u_(i-1); %velocity at start of panle
    [delx, thickhist] = ode45(@thickdash, [0, x(i)-x(i-1)], thick0);
    theta(i) = thickhist(end, 1);
    delta_e = thickhist(end, 2);
    Hes(i) = delta_e/theta(i);
    if Hes(i)>1.58 && ils ~= 0
        itr = i; %turbulent reattachment
    end
    if Hes(i) <1.46
        its = i; %turbulent separation
    end
end

if Hes(i)>=1.46
        H =(11*Hes(i) + 15)/(48*Hes(i)-59);
else
    H = 2.803;        
end 
if its ~= 0
    Hes(its+1:np) = Hes(i);
    theta(its+1:np) = theta(i).*(u_(i)./u_(i+1:np)).^(H+2);
end


disp(['re_l = ' num2str(re_l) ', dudx = ' num2str(dudx)])
if int ~= 0
    disp(['Natural transition at x = ' num2str(x(int)) ', Rethet/1000 = ' num2str(re_l * u_(int) * theta(int) / 1000)])
end
if ils ~= 0
    disp(['Separation at x = ' num2str(x(ils)) ', Rethet/1000 = ' num2str(re_l * u_(ils) * theta(ils) / 1000)])
end
if itr ~= 0
    disp(['Turbulent Reattachment at x = ' num2str(x(itr)) ', Rethet/1000 = ' num2str(re_l * u_(itr) * theta(itr) / 1000)])
end
if its ~= 0  
    disp(['Turbulent Seperation at x = ' num2str(x(its)) ', Rethet/1000 = ' num2str(re_l * u_(its) * theta(its) / 1000)])
end
if int==0 && ils ==0
    disp('No transition or separation')
end
disp(' ')

theta1 = theta;
Hes1 = Hes;

int1 =int;
ils1 = ils;
itr1 = itr;
its1 = its;

re_l = 1e6;

theta = zeros(1,np);  %theta squared at each position
Hes = zeros(1,np);
Hes(1) = 1.57258;


laminar = true;
integ = 0;
int = 0; %location of laminar->turbulent transition
ils = 0; %location of laminar seperation
itr = 0; %location of turbulent reattachment 
its = 0; %location of turbulent seperation
i=1;
u_ = linspace(1, 1+dudx, np);
while laminar && i<np
    i=i+1;
    integ = integ + ueintbit(x(i-1), u_(i-1), x(i), u_(i));
    theta(i) = sqrt((0.45/re_l) * u_(i)^(-6) * integ);
    m = -re_l * theta(i)^2  * dudx;
    H = thwaites_lookup(m);
    Hes(i) = laminar_He(H);
    re_thet = re_l * u_(i) * theta(i);
    if log(re_thet) >= 18.4*Hes(i) - 21.74
        laminar = false;    %transition
        int = i;

    elseif m>= 0.09
        laminar = false;
        ils = i;
        Hes(i) = 1.51509;
    end            
end
delta_e = Hes(i)*theta(i);

while its == 0 && i<np %loop for turbulent boundary layer
    i=i+1;
    thick0 = [theta(i-1), delta_e];
    u0 = u_(i-1); %velocity at start of panle
    [delx, thickhist] = ode45(@thickdash, [0, x(i)-x(i-1)], thick0);
    theta(i) = thickhist(end, 1);
    delta_e = thickhist(end, 2);
    Hes(i) = delta_e/theta(i);
    if Hes(i)>1.58 && ils ~= 0
        itr = i; %turbulent reattachment
    end
    if Hes(i) <1.46
        its = i; %turbulent separation
    end
end

if Hes(i)>=1.46
        H =(11*Hes(i) + 15)/(48*Hes(i)-59);
else
    H = 2.803;        
end 
if its ~= 0
    Hes(its+1:np) = Hes(i);
    theta(its+1:np) = theta(i).*(u_(i)./u_(i+1:np)).^(H+2);
end


disp(['re_l = ' num2str(re_l) ', dudx = ' num2str(dudx)])
if int ~= 0
    disp(['Natural transition at x = ' num2str(x(int)) ', Rethet/1000 = ' num2str(re_l * u_(int) * theta(int) / 1000)])
end
if ils ~= 0
    disp(['Separation at x = ' num2str(x(ils)) ', Rethet/1000 = ' num2str(re_l * u_(ils) * theta(ils) / 1000)])
end
if itr ~= 0
    disp(['Turbulent Reattachment at x = ' num2str(x(itr)) ', Rethet/1000 = ' num2str(re_l * u_(itr) * theta(itr) / 1000)])
end
if its ~= 0  
    disp(['Turbulent Seperation at x = ' num2str(x(its)) ', Rethet/1000 = ' num2str(re_l * u_(its) * theta(its) / 1000)])
end
if int==0 && ils ==0
    disp('No transition or separation')
end
disp(' ')

theta2 = theta;
Hes2 = Hes;


figure;
plot(x, theta1);
hold on;
plot(x, theta2);
hold off
legend('Re=1e6', 'Re=1e7');
xlabel('x/L');
ylabel('\theta/L');
if int ~= 0
    annotation('textarrow', [x(int)-0.3*x(int), x(int)],[theta1(int), theta1(int)] ,'String', 'Laminar Transition')
end
if ils ~= 0
    annotation('textarrow', [x(ils)-0.3*x(ils), x(ils)],[theta1(ils), theta1(ils)] ,'String', 'Laminar Separation')
end
if itr ~= 0
    annotation('textarrow', [x(itr)-0.3*x(itr), x(itr)],[theta1(itr), theta1(itr)] ,'String', 'Turbulent Reattachment')
end
if its ~= 0  
    annotation('textarrow', [x(its)-0.3*x(its), x(its)],[theta1(its), theta1(its)] ,'String', 'Turbulent Separation')
end
if int1 ~= 0
    annotation('textarrow', [x(int1)-0.3*x(int1), x(int1)],[theta1(int1), 0.2] ,'String', 'Laminar Transition')
end
if ils1 ~= 0
    annotation('textarrow', [x(ils1)-0.3*x(ils1), x(ils1)],[theta1(ils1), theta1(ils1)] ,'String', 'Laminar Separation')
end
if itr1 ~= 0
    annotation('textarrow', [x(itr1)-0.3*x(itr1), x(itr1)],[theta1(itr1), theta1(itr1)] ,'String', 'Turbulent Reattachment')
end
if its1 ~= 0  
    annotation('textarrow', [x(its1)-0.3*x(its1), x(its1)],[theta1(its1), theta1(its1)] ,'String', 'Turbulent Separation')
end



figure;
plot(x, Hes1);
hold on;
plot(x, Hes2);
hold off;
legend('Re=1e6', 'Re=1e7')
xlabel('x/L');
ylabel('He');
ylim([1.3,2.0]);
if int ~= 0
    xline(x(int), '-', 'Laminar Transition')
end
if ils ~= 0
    xline(x(ils), '-', 'Laminar Seperation')
end
if itr ~= 0
    xline(x(itr), '-', 'Turbulent Reattachment')
end
if its ~= 0  
    xline(x(its), '-', 'Turbulent Seperation')
end
