clear all;
close all;

my_labels = {'-30s', 'Pre10', '10s', '20s', '30s', '60s'};

datasets = {'SCD -> SC','SCD -> SC w/ antimycin'};
tableNames = {'pH5NoDex','pH5NoDexInhb'};
files = { 'ph5_no_dex.mat','ph5_no_dex_ihb.mat'};
colors = {[0 140/255 1],[222/255 140/255 0]};
linestyles={'-','-'};


load('annotationE221646');
annotationInfo = annotationE221646;
composite.names  = datasets;

%for the figure

%these are for the main figure (3)
ofinterest = [331 138 225 90 201 81 42 63 15 84 320 301 279 222 346];

%these are for the supplement S2; if you wish to see these, uncomment this
%line

%ofinterest = [129 84 15 98 80 63 81 37 88 139 42 110 40 78 341];

num_of_ions = length(ofinterest);
bins = [1:12:num_of_ions num_of_ions+1];
rownames = {};

for i=1:height(annotationInfo);
    rownames{i} = num2str(i);
end

annotationInfo.Properties.RowNames = rownames';

%load all of the data
for i=1:length(datasets)
    load(files{i});
    composite.data{i} = parsemyvar(myvar);
    A=table(cell(height(annotationInfo),1),'VariableNames',{tableNames{i}});
    A.Properties.RowNames = rownames';
    annotationInfo = join(annotationInfo,A,'Keys','RowNames');
end

foldChanges10 = zeros(1,num_of_ions);
foldChanges60 = zeros(1,num_of_ions);

bounceamount = zeros(1,num_of_ions);
drainamount = zeros(1,num_of_ions);

%declare indices of relevance

post_treatment = [3:6 9:12 15:18 21:24 27:30 33:36];
pre_treatment = [2 1 8 7 14 13 20 19 26 25 32 31];

for j=1:length(bins)-1
    figure;    
    for g=bins(j):bins(j+1)-1
        ion = ofinterest(g);
      
        subplot(3,4,g-bins(j)+1);
        hold on;
        maxVal = 0;
        minVal = 1000;
     
        ionInt = find(annotationInfo.ion  == ion);
        
        for k=1:length(datasets)
            composite.plotvect{k} = zeros(length(composite.data{k}(1).ODs),length(composite.data{k}));
            composite.rawvect{k} = composite.plotvect{k};
            for i=1:length(composite.data{k})
                composite.plotvect{k}(:,i) = log(composite.data{k}(i).data(:,ion) ./ composite.data{k}(1).ODs) * 1/log(2);
                composite.rawvect{k}(:,i) = composite.data{k}(i).data(:,ion);
            end    
            %to adjust by offset. comment out if undesired
            offset = mean(median(composite.plotvect{k}(:,1:2)));
            composite.plotvect{k}(:,:) = composite.plotvect{k}(:,:) -  mean(median(composite.plotvect{k}(:,1:2)));
            
            tran_pv = composite.plotvect{k}';
            
             %plot the SCD condition if the first data set
            if k==1
                firstSet=tran_pv(post_treatment(4:4:24));
            end
           
            
           plot(repmat([1:3 6],1,6),tran_pv(post_treatment),'.','Color',colors{k},'MarkerSize',2);
            hold on;
            

            [curve, goodness, output] = fit([2:6]',median(composite.plotvect{k}(:,2:end))','smoothingspline');
            errorbar([0:3 6],[0 median(composite.plotvect{k}(:,3:end))]', [0 std(composite.plotvect{k}(:,3:end))]' /sqrt(3),'Color',colors{k},'LineStyle',linestyles{k});
            if k==2
             plot(repmat([1 3],1,6),tran_pv(pre_treatment),'.','Color',[1 95/255 1],'MarkerSize',2);
        
             errorbar([0 1 3],[0 median(composite.plotvect{k}(:,[2 1]))]', [0 std(composite.plotvect{k}(:,[2 1]))]' /sqrt(3), 'Color', [1 95/255 1]);
             %calculate the significance
             secondSet=tran_pv(post_treatment(4:4:24));
             %calculating for significant changes
                [h, p] = ttest2(firstSet,secondSet,'Tail','both');
                %if it's a hit, plot it
                if h
                    text(6,max(ylim())-1,num2str(p));
                end
                
            end   
                legend off;
                %various metrics that we were keeping track of
                foldChanges10(ion) = mean(composite.plotvect{k}(:,3));
                foldChanges60(ion) = mean(composite.plotvect{k}(:,6));
                bounceamount(ion) = median(composite.plotvect{k}(:,5)) - median(composite.plotvect{k}(:,4));
                drainamount(ion) = median(composite.plotvect{k}(:,3));
                
                limits = ylim;
               

            maxVal = max([max(max(tran_pv)) maxVal]);
            minVal = min([min(min(tran_pv)) minVal]);
        end
        axis([0 7 floor(minVal) ceil(maxVal)]);
        ylabel('log2 change');
       
    
        if length(ionInt) < 3
            title([num2str(ion), ': ', annotationInfo.name{ionInt(1)}]);
        else
            title([num2str(ion), ': ', annotationInfo.name{ionInt(1)}]);
        end
        
        
    end
    %set(gcf,'PaperPositionMode','auto')
    %set(gcf, 'Position', get(0,'Screensize'));
    %print(strcat('adj_figure',num2str(j),'.png'),'-dpng','-r0');
end


