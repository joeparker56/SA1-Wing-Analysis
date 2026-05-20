function lhsmat = build_lhs_vec(xs, ys)
    np = length(xs)-1;    
    
    xa = xs(1:np)';
    xb = xs(2:np+1)';
    ya = ys(1:np)';
    yb = ys(2:np+1)';

    xc = xs(1:np);
    yc = ys(1:np);
    

    [infa, infb] = panelinf_vec(xa, ya, xb, yb, xc, yc);
    
    infa = infa';
    infb = infb';

    psip = [infa(:,1), infa(:,2:np) + infb(:,1:np-1), infb(:,np)];
    
    lhsmat = zeros(np+1, np+1);
    lhsmat(1:np-1,:) = psip(2:np,:) - psip(1:np-1,:);

    %set up the equation that gamma_np+1 = gamma_1 (circulation at last
    %point is the same as at the first point)
    lhsmat(np, 1) = 1;
    lhsmat(np, np+1) = -1;

    %set up the equation that gamma_1 = 0 (Kutta condition)
    lhsmat(np+1,1)=1;

end
