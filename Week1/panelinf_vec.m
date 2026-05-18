function [infa, infb] = panelinf_vec(xa, ya, xb, yb, x, y)
del = norm([xb-xa, yb-ya]);
tx = (xb-xa)/del;
ty = (yb-ya)/del;
nx_ = (yb-ya)/del;
ny_ = -(xb-xa)/del;

X = (x-xa)*tx + (y-ya)*ty;  %to find X for all x as x is a vector now, so we cant use the dot() function as we need to compute them element wise
Y = (x-xa)*nx_ + (y-ya)*ny_;

[infa, infb] = refpainf_vec(del, X, Y);
end