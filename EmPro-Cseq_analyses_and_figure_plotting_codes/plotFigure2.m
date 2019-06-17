% plot Figure 2 ABCDEF v3

close all
clear
clc
load expdata.mat
load sim1k.mat
load transcription_errors.mat

%%
expvaf = Exp.altreads./(Exp.depths + eps);
simvaf = Sim.altreads./(Sim.depths + eps);
simvafmean = mean(simvaf,2);
simvafstd = std(simvaf,0,2);
%%
pos = (1:9100)';

rdn25 = [922, 4318];
rdn58 = [4551, 4709];
rdn18 = [5070, 6870];
rdn5 = [8813, 8931];

% the loci to be plotted***
altreads = Exp.altreads;
ismut = expvaf > simvafmean + simvafstd*0 & ... 
    (...
    (pos >= rdn25(1) & pos <= rdn25(2)) | ...
    (pos >= rdn58(1) & pos <= rdn58(2)) | ...
    (pos >= rdn18(1) & pos <= rdn18(2)) | ...
    (pos >= rdn5(1) & pos <= rdn5(2)) ...
    );
altreads = altreads(ismut);
pv = pvalues(ismut);
nmut = nnz(ismut);

Ref = readtable('ref.txt', 'ReadVariableNames', false);
ref = Ref.Var1{1};
refmut = ref(ismut);

mutpos = find(ismut);
kscore = zeros(nmut,1);
neardist = zeros(nmut,1);

cc = zeros(nmut,3);
sz = ones(nmut,1)*30;
sz(pv < 0.05) = sz(pv < 0.05)*5;

for i = 1:nmut
    absdist = abs(mutpos - mutpos(i));
    absdist(i) = 99999;
    [neardist(i), nearix] = min(absdist);
    
    kscore(i) = neardist(i)/(altreads(i)*altreads(nearix));
    
    % mark color
    switch refmut(i)
        case 'A'
            cc(i,:) = hex2rgb('#CF3721');
        case 'T'
            cc(i,:) = hex2rgb('#F5BE41');
        case 'C'
            cc(i,:) = hex2rgb('#31A9B8');
        case 'G'
            cc(i,:) = hex2rgb('#258039');
        otherwise
            error('unknown nuc!')
    end
end
%%

zz = 24;

figure('Position', [0 0 2500 1250])

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Figure 2A: Histogram

subplot(3,4,1:2)
hold on

binnear = hist(neardist, 0:1:100);
bar(0:1:100, binnear, 0.618, 'FaceColor', [.7 .7 .7], 'EdgeColor', 'k')
xlim([-0.5 25 + 0.5])
xticks(0:1:50)

ylim([0 ceil(max(binnear)/10)*10])
xlabel('Distance to Nearest Mutant Position (Base)')
ylabel('Number of Mutant Positions')
set(gca,'tickdir','out','TickLength',[0.005 0.005],'fontsize',zz,'box','off')
hold off

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Figure 2C Kataegis

subplot(3,4,5:6)
hold on
alpha = 0.4;
refa = refmut == 'A';
scatter(mutpos(refa), kscore(refa), sz(refa,:), cc(refa,:), 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
reft = refmut == 'T';
scatter(mutpos(reft), kscore(reft), sz(reft,:), cc(reft,:), 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
refc = refmut == 'C';
scatter(mutpos(refc), kscore(refc), sz(refc,:), cc(refc,:), 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
refg = refmut == 'G';
scatter(mutpos(refg), kscore(refg), sz(refg,:), cc(refg,:), 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
legend({'A','T','C','G'},'location','southeast', 'Orientation', 'horizontal','fontsize',zz)

gy = 100;
dy = 80;
ty = 200;
gc = [.7 .7 .7];
plot1rectangle(rdn25(1),rdn25(2),gy,gy+dy,gc)
plot1rectangle(rdn58(1),rdn58(2),gy,gy+dy,gc)
plot1rectangle(rdn18(1),rdn18(2),gy,gy+dy,gc)
text(mean(rdn25), gy + ty, 'RDN25-2', 'FontSize',zz, 'HorizontalAlignment', 'center')
text(mean(rdn58), gy + ty, 'RDN58-2', 'FontSize',zz, 'HorizontalAlignment', 'center')
text(mean(rdn18), gy + ty, 'RDN18-2', 'FontSize',zz, 'HorizontalAlignment', 'center')

xticks([0 1000 2000 3000 4000 5000 6000 7000 8000 9000])
xticklabels({'460,000','461,000','462,000','463,000','464,000','465,000','466,000','467,000','468,000','469,000'})

xlabel('Chromosome XII')
ylabel('Kataegis Score')
xlim([500 7500])
ylim([1e-2 500])
set(gca,'tickdir','out','TickLength',[0.005 0.005],'fontsize',zz,'box','off', 'YMinorTick','off','YScale','log')
hold off


%%

load Fig2Cdata.mat

xb = 26;

E = expvaf';
S = simfreq';

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Figure 2B

subplot(3,4,9:10)
hold on
aboxplot(E(:,1:xb), 'colormap', hex2rgb('#ED5752'), 'outliermarkersize', 10)
aboxplot(S(:,1:xb), 'colormap', hex2rgb('#92AAC7'), 'outliermarkersize', 10, 'FaceAlpha', 0.5)

legend({'Experiment','Simulation'}, 'Location','north','Orientation','horizontal')

sigexppos = find(pvalues(1:xb) < 0.05);

xticks(1:xb)
xtbplain = arrayfun(@num2str, 0:xb-1, 'UniformOutput', false);
xtbcolor = xtbplain;
for i = 1:length(sigexppos)
    xtbcolor{i} = ['\color{red}', xtbplain{i}];
end

xticklabels(xtbcolor)

xlabel('Distance to 3-Prime End (Base)')
ylabel('Error Rate')
ylim([0 2.5e-4])
set(gca,'tickdir','out','TickLength',[0.005 0.005],'fontsize',zz,'box','off')
hold off

%%

C = readtable('Figure2D/ReadCount_3pr_1');
M = readtable('Figure2D/muta_count_at_1_mid.txt');

pos = 1:9100;
pos = pos';
n = length(pos);

dep3 = zeros(n,1);
mut3 = zeros(n,1);

for i = 1:size(C,1)
    dep3(C.position_in_genone(i)) = C.read_count_3pr(i);
end

for i = 1:size(M,1)
    mut3(M.position_in_genome(i)) = M.mutation_count(i);
end

rdn25 = [923, 4318];
rdn58 = [4551, 4708];
rdn18 = [5070, 6869];
rdn5 = [8813, 8931];
gy = -0.15;
dy = 0.3;
gc = [.7 .7 .7]; % color

logdep = log10(dep3);

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Figure 2D

subplot(3,4,[3,4,7,8])
hold on
bar(pos, logdep, 10, 'FaceColor', hex2rgb('#3A5199')); % blue
bar(pos, -mut3, 10, 'FaceColor', hex2rgb('#AF1C1C')); % red

plot1rectangle(rdn25(1),rdn25(2),gy,gy+dy,gc)
plot1rectangle(rdn58(1),rdn58(2),gy,gy+dy,gc)
plot1rectangle(rdn18(1),rdn18(2),gy,gy+dy,gc)

xlabel('Chromosome XII')
xlim([500 7500])
xticks(500:500:7500)
xticklabels({'','461,000','','462,000','','463,000','','464,000','','465,000','','466,000','','467,000',''})

% Bottom panel

T = readtable('Figure2D/muta_at_1_start_end_mutapos_bottom.txt','ReadVariableNames',false);
T.Properties.VariableNames = {'start_position', 'end_position', 'mutation_position'};

numtrs = size(T,1);

D = zeros(numtrs,n);

for i = 1:numtrs
    for k = 1:numtrs
        if sum(D(k,T.start_position(i):T.end_position(i))) == 0
            D(k,T.start_position(i):T.end_position(i)) = i;
            break
        end
    end
end

rmrows = sum(D,2) == 0;
D(rmrows,:) = [];

rr = zeros(numtrs,1);
for i = 1:size(D,1)
    for j = 1:size(D,2)
        if D(i,j) > 0 && rr(D(i,j)) == 0
            rr(D(i,j)) = i;
        end
    end
end

rng(1)
rr = rr + rand(size(rr));

for i = 1:numtrs
    line([T.start_position(i) T.end_position(i)],[-rr(i)/2-2 -rr(i)/2-2],...
        'color',hex2rgb('#3A5199'),'linewidth',5)
    mpos = strsplit(T.mutation_position{i}, ',');
    npos = -ones(1,length(mpos));
    for k = 1:length(mpos)
        npos(k) = str2double(mpos{k});
        if npos(k) == 0
            ss = 'o';
            cc = hex2rgb('#AF1C1C');
            sz = 6;
            plot(T.start_position(i) + npos(k), -rr(i)/2-2, ss,...
                'MarkerFaceColor',cc,'MarkerSize',sz,'markeredgecolor','none')
        end
    end
end

ylim([-5.5 4.5])
yticks(-2:4)
yticklabels({'2','1','0','10','100','1000','10000'})

set(gca,'tickdir','out','TickLength',[0.005 0.005],'fontsize',zz,'box','off')
hold off

%%

subplot(3,4,11)
hold on
bar(pos, logdep, 1, 'FaceColor', hex2rgb('#3A5199')); % blue
bar(pos, -mut3, 1, 'FaceColor', hex2rgb('#AF1C1C')); % red

plot1rectangle(rdn25(1),rdn25(2),gy,gy+dy,gc)
plot1rectangle(rdn58(1),rdn58(2),gy,gy+dy,gc)
plot1rectangle(rdn18(1),rdn18(2),gy,gy+dy,gc)

xlabel('Chromosome XII')
xlim([2760, 2790])
xticks([2760, 2770, 2780, 2790])
xticklabels({'462,760','462,770','462,780','462,790'})

% Bottom panel

numtrs = size(T,1);

D = zeros(numtrs,n);

for i = 1:numtrs
    for k = 1:numtrs
        if sum(D(k,T.start_position(i):T.end_position(i))) == 0
            D(k,T.start_position(i):T.end_position(i)) = i;
            break
        end
    end
end

rmrows = sum(D,2) == 0;
D(rmrows,:) = [];

rr = zeros(numtrs,1);
for i = 1:size(D,1)
    for j = 1:size(D,2)
        if D(i,j) > 0 && rr(D(i,j)) == 0
            rr(D(i,j)) = i;
        end
    end
end

rng(1)
rr = rr + rand(size(rr));

for i = 1:numtrs
    line([T.start_position(i) T.end_position(i)],[-rr(i)/2-2 -rr(i)/2-2],...
        'color',hex2rgb('#3A5199'),'linewidth',5)
    mpos = strsplit(T.mutation_position{i}, ',');
    npos = -ones(1,length(mpos));
    for k = 1:length(mpos)
        npos(k) = str2double(mpos{k});
        if npos(k) == 0
            ss = 'o';
            cc = hex2rgb('#AF1C1C');
            sz = 10;
            plot(T.start_position(i) + npos(k), -rr(i)/2-2, ss,...
                'MarkerFaceColor',cc,'MarkerSize',sz,'markeredgecolor','none')
        end
    end
end

ylim([-5.5 4.5])
yticks(-2:4)
yticklabels({'2','1','0','10','100','1000','10000'})

set(gca,'tickdir','out','TickLength',[0.01 0.01],'fontsize',zz,'box','off')
hold off

%%

subplot(3,4,12)
hold on
bar(pos, logdep, 1, 'FaceColor', hex2rgb('#3A5199')); % blue
bar(pos, -mut3, 1, 'FaceColor', hex2rgb('#AF1C1C')); % red

plot1rectangle(rdn25(1),rdn25(2),gy,gy+dy,gc)
plot1rectangle(rdn58(1),rdn58(2),gy,gy+dy,gc)
plot1rectangle(rdn18(1),rdn18(2),gy,gy+dy,gc)

xlabel('Chromosome XII')
xlim([5300, 5330])
xticks([5300, 5310, 5320, 5330])
xticklabels({'465,300','465,310','465,320','465,330'})

% Bottom panel

numtrs = size(T,1);

D = zeros(numtrs,n);

for i = 1:numtrs
    for k = 1:numtrs
        if sum(D(k,T.start_position(i):T.end_position(i))) == 0
            D(k,T.start_position(i):T.end_position(i)) = i;
            break
        end
    end
end

rmrows = sum(D,2) == 0;
D(rmrows,:) = [];

rr = zeros(numtrs,1);
for i = 1:size(D,1)
    for j = 1:size(D,2)
        if D(i,j) > 0 && rr(D(i,j)) == 0
            rr(D(i,j)) = i;
        end
    end
end

rng(1)
rr = rr + rand(size(rr));

for i = 1:numtrs
    line([T.start_position(i) T.end_position(i)],[-rr(i)/2-2 -rr(i)/2-2],...
        'color',hex2rgb('#3A5199'),'linewidth',5)
    mpos = strsplit(T.mutation_position{i}, ',');
    npos = -ones(1,length(mpos));
    for k = 1:length(mpos)
        npos(k) = str2double(mpos{k});
        if npos(k) == 0
            ss = 'o';
            cc = hex2rgb('#AF1C1C');
            sz = 10;
            plot(T.start_position(i) + npos(k), -rr(i)/2-2, ss,...
                'MarkerFaceColor',cc,'MarkerSize',sz,'markeredgecolor','none')
        end
    end
end

ylim([-5.5 4.5])
yticks(-2:4)
yticklabels({'2','1','0','10','100','1000','10000'})

set(gca,'tickdir','out','TickLength',[0.01 0.01],'fontsize',zz,'box','off')
hold off
