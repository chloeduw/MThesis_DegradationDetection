function deta = deta(a,b,c,x) %derivative of the electrode overpotentials
    num = -a*b*c*(1-2*x);
    denom = 2*c.*x.*(1-x).*sqrt(c^2.*x.*(1-x)+b^2);
    deta = num./denom;
end