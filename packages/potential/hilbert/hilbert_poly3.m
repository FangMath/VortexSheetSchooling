function pot = hilbert_poly3(a, w)
    % p(x) = a(1)x^n + a(2)x^{n-1} + ...
    % w: target points
    a = fliplr(a);

if length(a) ~= 4
    error('a should be length 4');
end

    lgd_coef = [a(1)+a(3)/3,...
                3*a(4)/5+a(2),...
                2*a(3)/3,...
                2*a(4)/5];

    y(1) = Q0(w);
    y(2) = Q1(w);
    y(3) = Q2(w);
    y(4) = Q3(w);
    pot = sum(y.*lgd_coef);

function y = Q0(z)
   y = complex(0,0);
   ztmp = (1+z)./(z-1);
   y = 0.5*log(ztmp);
end

function y = Q1(z)
   y= complex(0,0);
   ztmp = (1+z)./(z-1);
   y = 0.5*z*log(ztmp) - 1;
end

function y = Q2(z)
   y = complex(0,0);
   ztmp = (1+z)./(z-1);
   y = 0.25*(3*z.^2 -1)*log(ztmp) - 3*z./2;
end

function y = Q3(z)
   y = complex(0,0);
   ztmp = (1+z)./(z-1);
   y = 0.25*(5*z.^3 - 3*z)*log(ztmp) -2.5*z.^2 + 2/3;
end

end

