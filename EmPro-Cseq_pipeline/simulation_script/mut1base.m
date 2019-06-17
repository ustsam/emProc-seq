function bm = mut1base(bw)

N = {'A','T','G','C'};

if ~isempty(strcmp(N,bw))
    N(strcmp(N,bw)) = [];
else
    error('unknown input base!')
end

bm = N(randi(length(N),1));
bm = bm{1};