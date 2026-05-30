function lhsmat = build_lhs_vec(xs, ys)
    np = length(xs)-1;
    psip = zeros(np, np+1);
    infa = zeros(np, np+1);
    infb = zeros(np, np+1);

    
    [infa(1:np,1), infb(1:np,1)] = panelinf_vec(xs(1), ys(1), xs(2), ys(2), xs(1:np), ys(1:np));
    psip(1:np, 1) = infa(1:np,1);
    
    for j=2:np % main loop to generate psip matrix
        [infa(1:np, j), infb(1:np,j)] = panelinf_vec(xs(j), ys(j), xs(j+1), ys(j+1), xs(1:np), ys(1:np));
        psip(1:np,j) = infa(1:np,j) + infb(1:np,j-1);
    end

    %special case for j = np+1
    psip(1:np,np+1) = infb(1:np,j-1);
    lhsmat = zeros(np+1, np+1);
    lhsmat(1:np-1, :) = psip(2:np,:)-psip(1:np-1,:);

    %set up the equation that gamma_np+1 = gamma_1
    lhsmat(np, 1) = 1;
    lhsmat(np, np+1) = -1;

    %set up the equation that gamma_1 = 0
    lhsmat(np+1,1)=1;

end
