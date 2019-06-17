function plotregion_shaded(vaf0, vaf1, region, pv)

% region = 2000:3000;

position = 1:9100;
unplotted = setdiff(position,region);

vaf0 = vaf0'; % from pos*sim to sim*pos
vaf1 = vaf1';

vaf0reg = vaf0;
vaf0reg(:,unplotted) = 0;

vaf1reg = vaf1;
vaf1reg(:,unplotted) = 0;

pvsigreg = pv;
pvsigreg(unplotted) = 1;

shadedErrorBar(position,mean(vaf0reg),std(vaf0reg),'lineprops','b');
shadedErrorBar(position,mean(vaf1reg),std(vaf1reg),'lineprops','r');

position = position';

vaf1ub = mean(vaf1reg,1)' + std(vaf1reg)';

P = table(position, pvsigreg, vaf1ub);
P.Properties.VariableNames = {'pos','pval','vaf1ub'};
P0 = P(P.pval < 0.05,:);
plot(P0.pos, P0.vaf1ub, 'o','MarkerEdgeColor','k', 'MarkerFaceColor','r')
