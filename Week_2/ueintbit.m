function f = ueintbit(xa, ua, xb, ub)
    ubar = (ua+ub)/2;
    du = ub-ua;
    dx = xb-xa;
    f = (ubar^5 + (5/6)*ubar^3 * du^2 + (1/16) * ubar * du^4)*dx;
end

