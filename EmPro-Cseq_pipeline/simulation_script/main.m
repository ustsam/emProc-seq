% simulation of uniformly distributed errors

close all
clear
clc

%% Generate a random id
rng('shuffle')
symb1 = 'A':'Z';
symb2 = '0':'9';
nsmb1 = length(symb1);
nsmb2 = length(symb2);
id = [symb1(randperm(nsmb1,2)),symb2(randperm(nsmb2,4))];
fid = fopen(['sim',id,'.fastq.txt'], 'w');

%% Read reference rDNA unit
R = readtable('rDNA1.fasta.txt','ReadVariableNames',false);
R = R.Var1{2};
nr = length(R);

%% Start site and length
SL = readtable('startlengths.txt','ReadVariableNames',false);
SL.Properties.VariableNames = {'startsite','readlength'};

nt = size(SL,1);
s = SL.startsite;
e = SL.readlength;

%% Mutation rate
mutrate = 41881/66715376; % 0.0015

%% Generate transcripts with mutations and circularization

% Initialization
tr = cell(nt,1); % raw transcript, no circularization
trhead = cell(nt,1);

cirhead = cell(nt,1);
cirtr = cell(nt,1); % sequencing fragment with tandom repeats
cirqua = cell(nt,1);
reldist = zeros(nt,1); % relative distance of variant to 5' end
cr = 0;

aa = zeros(nt,1);

for i = 1:nt
    tr{i} = R(s(i):(s(i)+e(i)-1));
    trhead{i} = ['rDNA:',num2str(s(i)),'-',num2str(s(i)+e(i)-1)];
    
    % Mutating:
    pos = 1:e(i);
    mutix = rand(size(pos));
    [~,ai] = min(mutix);
    aa(i) = ai/e(i);
    nmut = nnz(mutix < mutrate);
    mutpos = pos(mutix < mutrate);
    muthead = [];
    if nmut > 0
        for k = 1:nmut
            cr = cr + 1;
            reldist(cr) = mutpos(k)/e(i);
            bw = tr{i}(mutpos(k));
            bm = mut1base(bw);
            tr{i}(mutpos(k)) = bm;
            muthead = strcat(muthead,'$mut:',num2str(s(i) + mutpos(k) - 1),':',bw,'>',bm);
        end
    end
    
    % Circularizing:
    cirpos = randi(e(i),1);
    cir1 = [tr{i}(cirpos:end),tr{i}(1:(cirpos-1))];
    cirtr{i} = [cir1,cir1,cir1(1:randi(e(i),1))]; % 2.x repeats
    cirhead{i} = strcat('@',trhead{i},'#cir:',num2str(s(i) + cirpos - 1),muthead);
    cirqua{i} = repmat('H',1,length(cirtr{i}));
    
    % Printing:
    fprintf(fid, '%s\n', cirhead{i});
    fprintf(fid, '%s\n', cirtr{i});
    fprintf(fid, '+\n');
    fprintf(fid, '%s\n', cirqua{i});
end
%%
reldist = reldist(1:cr);
%writetable(array2table(reldist), ['fastqs/',id,'.reldist.txt'],'WriteVariableNames',false)
%%
fclose(fid);