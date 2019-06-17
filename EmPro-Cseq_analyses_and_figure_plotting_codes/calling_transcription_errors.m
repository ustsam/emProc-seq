% transcriptional error calling at each position

clear
clc

load expdata.mat
load sim1k.mat

% Define exon regions:
rdn25 = [923, 4318];
rdn58 = [4551, 4708];
rdn18 = [5070, 6869];
rdn5 = [8813, 8931];

exonix = 1:nref;
exonix = exonix';
exonix(exonix >= rdn25(1) & exonix <= rdn25(2)) = 0;
exonix(exonix >= rdn58(1) & exonix <= rdn58(2)) = 0;
exonix(exonix >= rdn18(1) & exonix <= rdn18(2)) = 0;
exonix(exonix >= rdn5(1) & exonix <= rdn5(2)) = 0;

%% Call transcriptional errors

expvaf = Exp.altreads./(Exp.depths + eps);

calledpos = Exp.depths > 100 & exonix == 0;

binovaf1 = zeros(nref,nsim);

rng('shuffle')
for i = 1:nref
    if calledpos(i)
        disp(['calling ',num2str(i)])
        binovaf1(i,:) = binornd(Exp.depths(i), Exp.altreads(i)/Exp.depths(i), [1, nsim])/Exp.depths(i); 
    end
end

% Calculate p-value

simvaf = Sim.altreads./(Sim.depths + eps);

pvalues = ones(nref,1);

for i = 1:nref
    if calledpos(i)
        pvalues(i) = diffvaf(simvaf(i,:), binovaf1(i,:), 'paired');
    end
end

ix = (1:nref)';
meanVAFsim = mean(simvaf,2);

W = table(ix, Exp.altreads, Exp.depths, expvaf, meanVAFsim, pvalues);

%%
W.Properties.VariableNames = {'position','altered_reads','total_depth','VAF_exp','Mean_VAF_sim','p_value'};
W = sortrows(W,'p_value','ascend');
writetable(W,'Supp_Table1.txt','Delimiter','\t')

save('transcription_errors.mat','exonix','simvaf','binovaf1','pvalues','W')









