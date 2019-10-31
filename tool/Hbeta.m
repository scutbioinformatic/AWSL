function P = Hbeta(D, beta)
D = (D-min(D))/(max(D) - min(D)+eps);
P = exp(D * beta);
sumP = sum(P);
P = P / sumP;
end
