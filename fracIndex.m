function fracIndx = fracIndex(array,x)
fracIndx = zeros(1,length(x));
for idx = 1:length(x)
    if x(idx) >= array(length(array))
        fracIndx(idx) = length(array);
    elseif x(idx) <= array(1)
        fracIndx(idx) = 1;
    else
        a = find(array <= x(idx));
        a = a(length(a));
        b = find(array > x(idx));
        b = b(1);
        fracIndx(idx) = a+(x(idx)-array(a))/(array(b)-array(a));
    end
end

