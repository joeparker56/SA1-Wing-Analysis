function [int ils itr its delstar theta] = bl_solv_old(x,cp)
    global Re dudx u0
    np = length(x)+1;
    delstar = zeros(1,np);
    ue = [1,sqrt(1-cp)];
    xdash = [0,x];
    dudxs = (ue(2:np) - ue(1:np-1))./(xdash(2:np)-xdash(1:np-1));
    theta = zeros(1,np);
    Hes = zeros(1,np);
    Hes(1)=1.57258;
    laminar = true;
    integ=0;
    int=0;
    ils=0;
    itr=0;
    its=0;

    %length(dudxs);

    dudx = dudxs(1);
    integ = integ + ueintbit(xdash(1), ue(1), xdash(2), ue(2));
    theta(2) = sqrt((0.45/Re) * ue(2)^(-6) * integ);
    m = - Re * theta(2)^2  * dudx;
    H = thwaites_lookup(m);
    delstar(2) = H*theta(2);
    Hes(2) = laminar_He(H);
    i=2;

    while laminar && i<np
        i=i+1;
        dudx = dudxs(i-1);  %calculating theta i , so need to use dudx from the region before we get to i
        integ = integ + ueintbit(xdash(i-1), ue(i-1), xdash(i), ue(i));
        theta(i) = sqrt((0.45/Re) * ue(i)^(-6) * integ);
        m = -Re * theta(i)^2  * dudx;
        H = thwaites_lookup(m);
        delstar(i) = H*theta(i);
        Hes(i) = laminar_He(H);
        re_thet = Re * ue(i) * theta(i);

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
        u0 = ue(i-1); %velocity at start of panle
        dudx = dudxs(i-1);
        if i ~= np
            [delx, thickhist] = ode45(@thickdash, [0, xdash(i+1)-xdash(i)], thick0);
        else
            [delx, thickhist] = ode45(@thickdash, [0, xdash(i)-xdash(i-1)], thick0);
        end
        theta(i) = thickhist(end, 1);
        delta_e = thickhist(end, 2);
        Hes(i) = delta_e/theta(i);
        if Hes(i)>=1.46
            H =(11*Hes(i) + 15)/(48*Hes(i)-59);
        else
            H = 2.803;        
        end 
        delstar(i) = H*theta(i);        
        if Hes(i)>1.58 && ils ~= 0 &&itr==0
            itr = i; %turbulent reattachment
        end
        if Hes(i) <1.46
            its = i; %turbulent separation
        end
    end
    
    if its ~= 0
        Hes(its+1:np) = Hes(i);
        theta(its+1:np) = theta(i).*(ue(i)./ue(i+1:np)).^(H+2);
        for i=its+1:np
            i = its+1;
            if Hes(i)>=1.46
                H =(11*Hes(i) + 15)/(48*Hes(i)-59);
            else
                H = 2.803;        
            end 
            delstar(i) = H*theta(i);
        end

    end
    if int ~=0
        int = int-1;
    end
    if ils ~= 0      
        ils = ils-1;
    end
    if itr ~= 0 
        itr  = itr-1;
    end
    if its ~= 0 
        its = its-1;
    end
    theta = theta(2:end);
    delstar = delstar(2:end);
end
