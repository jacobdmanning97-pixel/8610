i0 = log(11/10);
for i = 1:20
    i0 = -10*i0 + 1/i;
end

fun = @(x) x.^20./(x+10);

i20 = integral(fun,0,1);

i0
i20