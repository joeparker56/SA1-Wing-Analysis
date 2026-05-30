function lhsmat = build_lhs(xs, ys)
    np = length(xs)-1;
    psip = zeros(np, np+1);
    infa = zeros(np, np+1);
    infb = zeros(np, np+1);

    for i=1:np %special case for j=1 
        [infa(i,1), infb(i,1)] = panelinf(xs(1), ys(1), xs(2), ys(2), xs(i), ys(i));
        psip(i, 1) = infa(i,1);
    end
    for j=2:np % main loop to generate psip matrix
        for i=1:np
            [infa(i, j), infb(i,j)] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
            psip(i,j) = infa(i,j) + infb(i,j-1);
        end
    end

    for i =1:np %special case for j = np+1
        psip(i,np+1) = infb(i,j-1);
    end
    lhsmat = zeros(np+1, np+1);
    for i=1:np-1
        for j=1:np+1
            lhsmat(i,j) = psip(i+1, j) - psip(i,j);
        end
    end
    %set up the equation that gamma_np+1 = gamma_1
    lhsmat(np, 1) = 1;
    lhsmat(np, np+1) = -1;

    %set up the equation that gamma_1 = 0
    lhsmat(np+1,1)=1;

end
