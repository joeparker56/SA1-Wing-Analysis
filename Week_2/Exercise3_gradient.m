clear;
close all;

np = 101;
laminar = true;
x = linspace(0,1,np); %length varies from 0->1

integ = 0; % value of the integral 
theta_2 = zeros(1,np);  %theta squared at each position
%re_ls=linspace(1e5,1e6,10000); searching for transition
dudx = -0.25;
found = false;
re_ls = [1e5, 1e4, 1e3];
for n=1:length(re_ls)
    if found
        break
    end
    
    re_l = re_ls(n);
    
    laminar = true;
    integ = 0;
    int = 0;
    ils = 0;
    i=1;
    u_ = linspace(1, 1+dudx, np);
    while laminar && i<np
        i=i+1;
        integ = integ + ueintbit(x(i-1), u_(i-1), x(i), u_(i));
        theta_2(i) = (0.45/re_l) * u_(i)^(-6) * integ;
        m = -re_l * theta_2(i)  * dudx;
        H = thwaites_lookup(m);
        He = laminar_He(H);
        re_thet = re_l * u_(i) * sqrt(theta_2(i));
        if log(re_thet) >= 18.4*He - 21.74
            laminar = false;    %transition
            int = i;

        elseif m>= 0.09
            laminar = false;
            ils = i;
        end            
    end
    if int ~= 0
        found =true;
        disp(['re_l = ' num2str(re_l) ', dudx = ' num2str(dudx)])
        %disp([num2str(re_ls(n-1)) ', ' num2str(re_l) ', ' num2str(re_l) ', ' num2str(re_ls(n+1))])        %surrounding reynolds numbers to check we have enough sig figs
        disp(['Natural transition at x = ' num2str(x(int)) ', Rethet/1000 = ' num2str(re_l * u_(int) * sqrt(theta_2(int)) / 1000)])
    elseif ils ~= 0
        disp(['re_l = ' num2str(re_l) ', dudx = ' num2str(dudx)])
        %disp([num2str(re_ls(n-1)) ', ' num2str(re_l) ', ' num2str(re_l) ', ' num2str(re_ls(n+1))])        

        disp(['Separation at x = ' num2str(x(ils)) ', Rethet/1000 = ' num2str(re_l * u_(ils) * sqrt(theta_2(ils)) / 1000)])
    else
        disp(['re_l = ' num2str(re_l) ', dudx = ' num2str(dudx)])
        disp('No transition or separation')
    end
    disp(' ')
    
end

%{
Exercise 2 Results:
re_l
     5000000

dudx
   -0.1000

x, re_thet/1000
    0.4800    1.0802

re_l
     5000000

dudx
     0

x, re_thet/1000
    0.6600    1.2870

re_l
    10000000

dudx
   -0.1000

x, re_thet/1000
    0.2900    1.1686

re_l
    10000000

dudx
     0

x, re_thet/1000
    0.3500    1.2902

re_l
    10000000

dudx
    0.1000

x, re_thet/1000
    0.5500    1.6452

re_l
    20000000

dudx
   -0.1000

x, re_thet/1000
    0.1700    1.2532

re_l
    20000000

dudx
     0

x, re_thet/1000
    0.1800    1.2906

re_l
    20000000

dudx
    0.1000

x, re_thet/1000
    0.2100    1.3973
%}

%{
Exercise 3 Results
re_l
      100000

dudx
   -0.2500

x, re_thet/1000
laminar seperation
    0.6400    0.1789

Seperation at 0.64 with Rethet 178.9051
   
re_l
       10000

dudx
   -0.2500

x, re_thet/1000
laminar seperation
    0.6400    0.0566

Seperation at 0.64 with Rethet 56.5748
   
re_l
        1000

dudx
   -0.2500

x, re_thet/1000
laminar seperation
    0.6400    0.0179

Seperation at 0.64 with Rethet 17.8905
%}