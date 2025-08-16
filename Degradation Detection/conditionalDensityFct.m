function f = conditionalDensityFct(NT,innov,S)
    f = zeros(1,NT);
    for i=1:NT
        f(i) = 1/((2*pi)^(1/2)*sqrt(det(S(i))))*exp(-1/2*innov(i)'*(S(i)\innov(i)));
    end
end