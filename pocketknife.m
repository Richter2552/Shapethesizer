
clear all
s=24;% define the number of shapes
n=1000;% define the number of iteration
rule = {'relaxed', 'strict'};
sym=[1,1;1,0;0,1;0,0];
% create a set of "s" number of shapes randomising "rules" and "symmetries"
for jj=1:s
genes(:,:,jj)=shape_proc(jj,[4,4],rule(randi(2)),sym(randi(4),:),0,0);
end
vec=1:s;
comb= nchoosek(vec,2);% % create all possible combinations of two from genes matrices 


for l=1:n
    %create a matrix indicating row numbers in which each element is
    %present in comb e.g., second element in genes is present in row numbers
    %row(2,:), third elements in row(3,:) in comb and so on
    for j=1:s
        row(j,:) = vertcat(find(ismember(comb(:,1),vec(j))),find(ismember(comb(:,2),vec(j))));
    end
    
    
    % take each pair of gene elements and evaluate similarity coefficient between them
    % with or without flipping and rotating them
    for i=1:size(comb,1)
        a= genes(:,:,comb(i,1));% first element
        b= genes(:,:,comb(i,2));% second element
        
       
        %evaluate similarity coefficient between two genes
        dv=(a-b); 
        dv(dv~=0)=1; % insert 1 in places where elements differ
        dv(:)=~dv; % replaces 1 and 0
        v(1)=sum(sum(dv))/16; % 
        
        %rotate one first gene counterclockwise 90 degrees and evaluate similarity 
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
    
   
    for ii=1:size(row,1)
        c(ii)=mean(vmax(1,row(ii,:)));
    end
    
    [cul,icul]=max(c);% identify the culptit and its index. Incase there are multiple it identifies the first culprit
    genes(:,:,icul)=shape_proc(1,[4,4],rule(randi(2)),sym(randi(4),:),0,0);% eliminate the culprit from the list
    
    
    clear row
    clear j
    clear ii
    clear c
end
% plot the iteration curve
x=1:n;
figure
plot(x,mavg,'b--o')
title('Change in the similarity coefficient among shapes with iterated replacement of members')
xlabel('No. of Iteration')
ylabel('Average Gene Value')
pl = gca;
pl.FontSize = 18;
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));

