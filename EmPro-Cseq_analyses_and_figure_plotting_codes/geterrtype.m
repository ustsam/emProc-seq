function errtypes = geterrtype(atypes, refseq, altype)
% count error types

errtypes = zeros(1,12);

n = length(refseq);

for i = 1:n
    if ~isempty(altype{i})        
        altnuclist = strsplit(altype{i}, '/');
        numalts = length(altnuclist);
        for k = 1:numalts
            if ~isempty(altnuclist{k})
                a1 = [refseq{i},'>',altnuclist{k}];
                errtypes(strcmp(atypes, a1)) = errtypes(strcmp(atypes, a1)) + 1;
            end
        end
    end
end
