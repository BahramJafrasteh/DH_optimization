clc
clear
close
file = 'Final_Results_EI_0.25.txt' ;
outf = 'AdBH_Esfordi_EI_0.25.txt';
load 'polygon_esfordi'
min_length_dh = 65;
max_length_dh = 200;
dist_bet_dh = 50;
ind_pi = 9;%9
ind_ei = 8;
medYh = 6;
medYs = 7;
ind_l = 10;


FinalResultsEsfordi = dlmread(file,' ',2,0);
FinalResultsEsfordi = FinalResultsEsfordi(:,2:end);

%% In polygon
[in,on] = inpolygon(FinalResultsEsfordi(:,1),FinalResultsEsfordi(:,2),polygon(:,1),polygon(:,2));

FinalResultsEsfordi = FinalResultsEsfordi(in,:);
%%



Collars = dlmread('Collars_esfordi.csv',',',0,0);
[a,ind] = sort(FinalResultsEsfordi(:,ind_pi));
sorted = FinalResultsEsfordi(ind,:);



ind = sorted(:,ind_pi) > 10000;
finalc = sorted(~ind,:);

[~,ind] = sort(-finalc(:,ind_pi));
finalc = finalc(ind,:);


for i = 1: size(finalc,1)
    dist = pdist2(Collars(:,1:2), finalc(i,1:2));
    dis(i) = min(dist);
end
dis = dis';
finalc = [finalc,dis];

ind_dis = finalc(:,end)>dist_bet_dh;
finalc_dis = finalc(ind_dis,:);
ind_length = finalc_dis(:,ind_l)<=max_length_dh & finalc_dis(:,ind_l)>min_length_dh;
finalc_length = finalc_dis(ind_length,:);


%dlmwrite('AdBH_Esfordi.txt', [finalc(1:10,:), (1:10)']);
M = [finalc_length(1:20,:), (1:20)'];
fid = fopen(outf,'wt');
fprintf(fid, '%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t  %s\t  %s\t  %s\t %s\n', 'x','y','z', 'az', 'dip', 'medYh', 'medYs', 'PI', 'EI', 'length', 'Dist_collar', 'rank');
% have a matrix M(N,3) ready to go
dlmwrite(outf, M,'delimiter', '\t', '-append')
fclose(fid);

