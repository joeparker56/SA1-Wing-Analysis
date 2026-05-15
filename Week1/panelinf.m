function [infa, infb] = panelinf(xa, ya, xb, yb, x, y)
    del = norm([(xb-xa), (yb-ya)]);
    r = [(x-xa), (y-ya)];
    t = [(xb-xa), (yb-ya)]/del;
    n = [(yb-ya), -(xb-xa)]/del;
    X = dot(r, t);
    Y = dot(r, n);
    [infa, infb] = refpainf(del, X, Y);
end 
