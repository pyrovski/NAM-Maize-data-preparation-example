%{
x must be sorted in increasing order
find index s.t. x(1:index) < el
if x >= el, return -1
%}

function [index] = binarySearch(x, el)
    high = length(x);
    low = 1;
    if el > x(high)
        index = high;
        return;
    elseif el <= x(low)
        index = -1;
        return;
    end
    while(low < high - 1)
        mid = round((low+high)/2);
        if(el <= x(mid))
            high = mid;
        else
            low = mid;
        end
    end
    if(el > x(low))
        index = low;
    else
        index = high;
    end
end
