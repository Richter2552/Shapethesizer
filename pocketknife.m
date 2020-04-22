clear all
s=24;% define the number of shapes
n=200;% define the number of iteration
rule = {'relaxed', 'strict'};
sym=[1,1;1,0;0,1;0,0];
%create a set of "s" number of shapes randomising "rules" and "symmetries"
for i=1:s
genes(:,:,i)=shape_proc(i,[4,4],rule(randi(2)),sym(randi(4),:),0,0);
end
vec=1:s;
comb= nchoosek(vec,2);% % create all possible combinations of two from genes matrices 


for it=1:n
    %create a matrix indicating row numbers in which each element is
    %present in comb e.g., second element in genes is present in row numbers
    %row(2,:), third elements in row(3,:) in comb and so on
    for j=1:s
        row(j,:) = vertcat(find(ismember(comb(:,1),vec(j))),find(ismember(comb(:,2),vec(j))));
    end
    
    
    % take each pair of gene elements and evaluate similarity coefficient between them
    % with or without flipping and rotating them
    for k=1:size(comb,1)
        a= genes(:,:,comb(k,1));% first element
        b= genes(:,:,comb(k,2));% second element
        
        % while calculating the similarity index among genes the microscopic features
        % (i.e., identity of individual primitive shapes) is ignored while the macrosopic
        % features (i.e., orientation of those primitive shapes) are
        % retained
        a(a == 20 | a == 30 | a == 40) = 0;
        a(a == 21 | a == 31 | a == 41) = 1;
        a(a == 22 | a == 32 | a == 42) = 2;
        a(a == 23 | a == 33 | a == 43) = 3;
        
        b(b == 20 | b == 30 | b == 40) = 0;
        b(b == 21 | b == 31 | b == 41) = 1;
        b(b == 22 | b == 32 | b == 42) = 2;
        b(b == 23 | b == 33 | b == 43) = 3;
        
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
        vmax(k)=max(v);
        
    end
    
    % calculate the average similarity coefficient among all the members of
    % the existing list
    mavg(it)=mean(vmax);
    mavg_temp_it=mean(vmax);
    mavg_temp=10;
    genes_temp_it=genes;
    
   %calculate the average contribution of each member of the existing list
   %in the similarity coeeficient
      for l=1:size(row,1)
         c(l)=mean(vmax(1,row(l,:)));
      end
    
    [cul,icul]=max(c);% identify the culptit and its index. Incase there are multiple it identifies the first culprit
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   while mavg_temp > mavg_temp_it
 
   genes_temp_it(:,:,icul)= shape_proc(1,[4,4],rule(randi(2)),sym(randi(4),:),0,0);
    
    
    for kk=1:size(comb,1)
        a= genes_temp_it(:,:,comb(kk,1));% first element
        b= genes_temp_it(:,:,comb(kk,2));% second element
        
        % while calculating the similarity index among genes the microscopic features
        % (i.e., identity of individual primitive shapes) is ignored while the macrosopic
        % features (i.e., orientation of those primitive shapes) are
        % retained
        a(a == 20 | a == 30 | a == 40) = 0;
        a(a == 21 | a == 31 | a == 41) = 1;
        a(a == 22 | a == 32 | a == 42) = 2;
        a(a == 23 | a == 33 | a == 43) = 3;
        
        b(b == 20 | b == 30 | b == 40) = 0;
        b(b == 21 | b == 31 | b == 41) = 1;
        b(b == 22 | b == 32 | b == 42) = 2;
        b(b == 23 | b == 33 | b == 43) = 3;
        
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
        
%         if v(1)==0 
%            continue
%         end  
        
        %calculate the maximum similarity coefficient with or without
        %rotating/flipping
        vmax_temp(kk)=max(v);
        clear v
    end
    
    % calculate the average similarity coefficient among all the members of
    % the existing list
    mavg_temp=mean(vmax_temp);
    
    %genes_temp_it(:,:,idul)=shape_proc(1,[4,4],rule(randi(2)),sym(randi(4),:),0,0);% replace the culprit from the list wit a randomly generated shape
    
   %clear row
    clear ii
    clear jj
    clear kk
    clear cc
   end
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   genes=genes_temp_it;
 
    clear row
    clear i
    clear j
    clear k
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

