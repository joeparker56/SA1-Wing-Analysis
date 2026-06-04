function [cl cd] = forces(circ,cp,delstarl,thetal,delstaru,thetau)
    cl = -2*circ;

    theta_te=thetal(length(thetal))+thetau(length(thetau));
    delstar=delstarl(length(delstarl))+delstaru(length(delstaru));
    He = delstar/theta_te;
    ue = sqrt(1-cp(end));
    theta_inf = theta_te*(ue)^((He+5)/2);

    cd = 2*theta_inf;
