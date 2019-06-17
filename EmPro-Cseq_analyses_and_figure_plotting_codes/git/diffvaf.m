function pv = diffvaf(v0,v1,comptype)

N = length(v0);
if length(v1) ~= N
    error('The input v0 and v1 must be in the same length!')
end

if strcmp(comptype,'direct')
    pv = nnz(v0 >= v1)/N;
elseif strcmp(comptype,'paired')
	s = 0;
    w = 0;
    pv = 1;
    for i = 1:N
        for j = 1:N
            w = w + 1;
            if v0(i) >= v1(j) % H0 is true
                s = s + 1;
            end
        end
    end

    if w > 0
        pv = s/w;
    end
    
else
    error('unknown p-value comparing method!')
end


