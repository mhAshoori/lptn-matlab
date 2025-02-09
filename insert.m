function x = insert(a , x , n, align)
%INSERT insert a scalar inside a vector
%   Detailed explanation goes here


if n == 1
    if align == +1
        x = [x(1);a;x(2:end)];
    else
        x = [a;x];

    end

elseif n==length(x)

    if align == +1
        x = [x;a];

    else
        x = [x(1:end-2);a;x(end-1:end)];
    end

else
    if align == +1
        x = [x(1:n); a; x(n+1:end)];
    else
        x = [x(1:n-1); a; x(n:end)];
    end

end

