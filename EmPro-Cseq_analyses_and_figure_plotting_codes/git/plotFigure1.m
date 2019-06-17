% Plot Figure 1

clc
clear
close all

load expdata.mat
load sim1k.mat

%% Figure 1B: PROseq error (3 batches) vs emProc error (5 batches)

proseq_error = [26303;27390;6746];
proseq_coverage = [74133277;64944918;18772664];

proseq_errormean = mean(proseq_error./proseq_coverage);
proseq_errorstd = std(proseq_error./proseq_coverage);

emproc_error = sum(L.altreads,1);
emproc_coverage = sum(L.depths,1);

emproc_errormean = mean(emproc_error./emproc_coverage);
emproc_errorstd = std(emproc_error./emproc_coverage);

%% Ploting
figure('Position',[0 0 1800 1200])
subplot(3,4,1)
hold on

bh1 = bar(1,proseq_errormean,0.4);
bh2 = bar(2,emproc_errormean,0.4);
bh1.FaceColor = hex2rgb('#E4E3DB');
bh2.FaceColor = hex2rgb('#E4E3DB');
xticks(1:2)
xticklabels({'Pro-seq','rPOLE-seq'})
ylabel('Mutational Frequency')

ber = errorbar([1,2],[proseq_errormean,emproc_errormean],[proseq_errorstd,emproc_errorstd], 'CapSize', 12, 'LineWidth', 1);
ber.Color = [0 0 0];
ber.LineStyle = 'none';
xlim([0.5 2.5])
set(gca,'tickdir','out','TickLength',[0.005 0.005],'fontsize',24,'box','off')
hold off

% Figure 1C

atypesDNA = {'T>A';...
            'T>G';...
            'T>C';...
            'A>T';...
            'A>G';...
            'A>C';...
            'G>T';...
            'G>A';...
            'G>C';...
            'C>T';...
            'C>A';...
            'C>G'};
        
atypesRNA = {'A>U';...
            'A>C';...
            'A>G';...
            'U>A';...
            'U>C';...
            'U>G';...
            'C>A';...
            'C>U';...
            'C>G';...
            'G>A';...
            'G>U';...
            'G>C'};
        
errtypes1 = geterrtype(atypesDNA, ref, L.altypes(:,1))/sum(L.depths(:,1));
errtypes2 = geterrtype(atypesDNA, ref, L.altypes(:,2))/sum(L.depths(:,2));
errtypes3 = geterrtype(atypesDNA, ref, L.altypes(:,3))/sum(L.depths(:,3));
errtypes4 = geterrtype(atypesDNA, ref, L.altypes(:,4))/sum(L.depths(:,4));
errtypes5 = geterrtype(atypesDNA, ref, L.altypes(:,5))/sum(L.depths(:,5));

errtypes = mean([errtypes1;errtypes2;errtypes3;errtypes4;errtypes5]);
stderrtypes = std([errtypes1;errtypes2;errtypes3;errtypes4;errtypes5]);

subplot(3,4,2:4)
hold on

bh1 = bar([1,2,3],errtypes([1,2,3]),0.6);
bh2 = bar([4,5,6],errtypes([4,5,6]),0.6);
bh3 = bar([7,8,9],errtypes([7,8,9]),0.6);
bh4 = bar([10,11,12],errtypes([10,11,12]),0.6);

bh1.FaceColor = hex2rgb('#FFBEBD');
bh2.FaceColor = hex2rgb('#EB8A3E');
bh3.FaceColor = hex2rgb('#EBB582');
bh4.FaceColor = hex2rgb('#ED8C72');

ber = errorbar(1:12,errtypes,stderrtypes, 'CapSize', 12, 'LineWidth', 1);
ber.Color = [0 0 0];
ber.LineStyle = 'none';

xticks(1:12)
xticklabels(atypesRNA)
%xtickangle(45)
xlim([0.4 12.6])
ylim([0 1.4e-5])
ylabel('Mutational Frequency')
set(gca,'tickdir','out','TickLength',[0.005 0.005],'fontsize',24,'box','off')
hold off

% Figure 1D: Landscape

subplot(3,4,5:8)
hold on

bar(1:nref,Exp.altreads, 10, 'EdgeColor', 'none', 'FaceColor', hex2rgb('#AF1C1C')); % red
bar(1:nref,-Exp.depths/5000, 10, 'EdgeColor', 'none', 'FaceColor', hex2rgb('#3A5199')); % blue

% rectangles for gene structure annotation:
rdn25 = [923, 4318];
rdn58 = [4551, 4708];
rdn18 = [5070, 6869];
rdn5 = [8813, 8931];
gy = -20;
gc = [.7 .7 .7];
plot1rectangle(rdn25(1),rdn25(2),gy,gy+2,gc)
plot1rectangle(rdn58(1),rdn58(2),gy,gy+2,gc)
plot1rectangle(rdn18(1),rdn18(2),gy,gy+2,gc)
plot1rectangle(rdn5(1),rdn5(2),gy,gy+2,gc)
fz = 24;
text(mean(rdn25), gy + 4, 'RDN25-2', 'FontSize',fz, 'HorizontalAlignment', 'center')
text(mean(rdn58), gy + 4, 'RDN58-2', 'FontSize',fz, 'HorizontalAlignment', 'center')
text(mean(rdn18), gy + 4, 'RDN18-2', 'FontSize',fz, 'HorizontalAlignment', 'center')
text(mean(rdn5), gy + 4, 'RDN5-2', 'FontSize',fz, 'HorizontalAlignment', 'center')

%
xlim([0 nref+1])
xticks(0:500:9000)
xticklabels({'460000','','461000','','462000','','463000','','464000','','465000','','466000','','467000','','468000','','469000'})

ylim([-21 11])
yticks(-20:5:10)
yticklabels({'100000','75000','50000','25000','0','5','10'})

xlabel('Chromosome XII')
ylabel('     Coverage       Errors      ')
xlim([0 nref+1])
set(gca,'tickdir','out','TickLength',[0.005 0.005],'fontsize',24,'box','off')
hold off

%% Figure E and F

load transcription_errors.mat
simvaf(exonix > 0,:) = 0;
binovaf1(exonix > 0,:) = 0;

pos = 1:9100;
depos = pos(Exp.depths > 100);

rdn25 = 923:4318;
rdn58 = 4551:4708;
rdn18 = 5070:6869;
rdn5 = 8813:8931;

yroof = 0.01;

subplot(3,4,9:12)
hold on
region = 800:7200;

r25 = intersect(region,rdn25);
r58 = intersect(region,rdn58);
r18 = intersect(region,rdn18);
r5 = intersect(region,rdn5);

u1 = union(r25, r58);
u2 = union(u1, r18);
u3 = union(u2, r5);

regionplotted = intersect(u3, depos);

plotregion_shaded(simvaf, binovaf1, regionplotted, pvalues);
xlim([min(region)-1, max(region)+1])
ylim([0 yroof])
xlabel('Chromosome XII')
ylabel('Mutational Frequency')
xt = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', strcat('46',xt))
set(gca,'tickdir','out','TickLength',[0.005 0.005],'fontsize',24,'box','off')
hold off

