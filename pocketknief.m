clear all
load sham_genes.mat % load a sham list of gene 


m=transpose(who);
n=17;% define number of iteration

for l=1:n
    
    comb= nchoosek(m,2);% create all possible combination of the matrices
    
    for j=1:size(m,2)
        row(j,:) = vertcat(find(ismember(comb(:,1),m(j))),find(ismember(comb(:,2),m(j))));
    end
    
    
    
    for i=1:size(comb,1)
        a=eval(comb{i,1});
        b=eval(comb{i,2});
        
        %evaluate similarity coefficient between two shapes
        dv=(a-b);
        dv(dv~=0)=1;
        dv(:)=~dv;
        v(1)=sum(sum(dv))/16;
        
        %rotate one image counterclockwise 90 degrees and evaluate similarity 
        a=rot90(a,1);
        dv=(a-b);
        dv(dv~=0)=1;
        dv(:)=~dv;
        v(2)=sum(sum(dv))/16;
        
        %rotate one image counterclockwise 180 degrees and evaluate similarity 
        a=rot90(a,2);
        dv=(a-b);
        dv(dv~=0)=1;
        dv(:)=~dv;
        v(3)=sum(sum(dv))/16;
        
         %rotate one image counterclockwise 270 degrees and evaluate similarity 
        a=rot90(a,3);
        dv=(a-b);
        dv(dv~=0)=1;
        dv(:)=~dv;
        v(4)=sum(sum(dv))/16;
        
         %flip one image horizontaly and evaluate similarity 
        a=fliplr(a);
        dv=(a-b);
        dv(dv~=0)=1;
        dv(:)=~dv;
        v(5)=sum(sum(dv))/16;
        
        %flip one image vertically and evaluate similarity
        a=flipud(a);
        dv=(a-b);
        dv(dv~=0)=1;
        dv(:)=~dv;
        v(6)=sum(sum(dv))/16;
        
        %calculate the maximum similarity coefficient with or without
        %rotating/flipping
        vmax(i)=max(v);
        
    end
    
    % calculate the average similarity coefficient among all the members of
    % the existing list
    mavg(l)=mean(vmax);
    
    x=[mean(vmax(1,row(1,:))),mean(vmax(1,row(2,:))),mean(vmax(1,row(3,:))),mean(vmax(1,row(4,:)))];
    
    
    [cul,icul]=max(x);% identify the culptit and its index. Incase there are multiple it identifies the first culprit
    m(icul)=[];% eliminate the culprit from the list
    
    clear row
    clear j
end
% plot the iteration curve
x=1:n;
figure
plot(x,mavg,'b--o')
title('Decrease in the similarity among shapes with iterated removal of members')
xlabel('No. of Iteration')
ylabel('Average Gene Value')
pl = gca;
pl.FontSize = 18;
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
