function [gene] = shape_proc(num,sz,rule,sym,draw,prnt)

% Procedural random shape generator

if any(sym & rem(sz,2))
    error('for this symmetry the size on this axis must be equal');
end

% Define primitives

sqr.coord = [0 0 1 1; 0 1 1 0];
sqr.area = 1;
sqr.lego = [1 1 1 1]; % specifies 'build-ability' of left,top,right,bottom

tri.coord = [0 0 1 0;0 1 0 0];
tri.area = 0.5;
tri.lego = [1 0 0 1];

cv.coord = [cos(linspace(0, pi/2)) 0; sin(linspace(0, pi/2)) 0];
cv.area = 1 - pi/4;
cv.lego = [1 0 0 1];

cc.coord = [1-cos(linspace(0, pi/2)) 0; 1-sin(linspace(0, pi/2)) 0];
cc.area = pi/4;
cc.lego = [1 0 0 1];

shape(1) = sqr;
shape(2) = tri;
shape(3) = cv;
shape(4) = cc;

%% start procedure

target = [0 0]; % initial shape origin (bottom/lefthand)
ang = [0 90 180 270];

ii = 1;
shp = randsample(1:4,1);
rot_ang = randsample(ang,1);

cell_list(ii,:) = target;
shape_list(ii) = rot(shape(shp),shp,rot_ang);
gface{ii} = find(shape_list(ii).lego);
lego(ii,:) = shape_list(ii).lego;

c3 = 0;
for kk=gface{ii}
    c3 = c3 + 1;
    if kk == 1
        target_buff(c3,:) = [target(1)-1 target(2)];
    elseif kk == 2
        target_buff(c3,:) = [target(1) target(2)+1];
    elseif kk == 3
        target_buff(c3,:) = [target(1)+1 target(2)];
    elseif kk == 4
        target_buff(c3,:) = [target(1) target(2)-1];
    end
end

% enforce symmetry by adjusting bounding box
sz_r = sz;
if sym(1) && ~sym(2) % horizontal symmetry
    sz_r(2) = sz_r(2)/2;
elseif sym(2) && ~sym(1) % vertical symmetry
    sz_r(1) = sz_r(1)/2;
elseif sym(1) && sym(2) % dual symmetry
    sz_r(1) = sz_r(1)/2;
    sz_r(2) = sz_r(2)/2;
end

% close all
% figure
% patch(shape_list(ii).coord(1,:) + target(1),shape_list(ii).coord(2,:) + target(2),'k')

% Main drawing loop after initialization
while ~isempty(target_buff)
    target_buff = target_buff(randperm(size(target_buff,1)),:); % randomize order for even growth
    target_buff(ismember(target_buff,cell_list,'rows'),:) = []; % remove full possible targets
    % enforce symmetry condition
    if sym(1) && ~sym(2) % horizontal symmetry
        target_buff(target_buff(:,2)<0,:) = []; % remove targets across h-axis
    elseif sym(2) && ~sym(1) % vertical symmetry
        target_buff(target_buff(:,1)<0,:) = []; % remove targets across v-axis
    elseif sym(1) && sym(2) % dual symmetry
        target_buff(target_buff(:,1)<0,:) = []; % remove targets across v-axis
        target_buff(target_buff(:,2)<0,:) = []; % remove targets across h-axis
    end
    if ~isempty(target_buff)
        target = target_buff(1,:); % choose next position
        % check that max size is not exceeded
        while (max([cell_list(:,1); target(1)]) - min([cell_list(:,1); target(1)])) == sz_r(1) || (max([cell_list(:,2); target(2)]) - min([cell_list(:,2); target(2)])) == sz_r(2)
            target_buff(1,:) = [];
            if isempty(target_buff)
                gene = shape_proc(num,sz,rule,sym,draw,prnt)
                return
            end
            target = target_buff(1,:); % choose next position
        end
        
        % query compatibility conditions of new target cell (so conditions
        % outside)
        t_coord = [target(1) - 1 target(2) + 0; target(1) + 0 target(2) + 1; target(1) + 1 target(2) + 0; target(1) + 0 target(2) - 1];
        t_filled = ismember(t_coord,cell_list,'rows')'; % find filled cells (L/T/R/B)
        
        lego_inv = lego(:,[3 4 1 2]); % Needs to be re-ordered to R/B/L/T
        t_comp = ones(1,4);
        for jj = 1:4
            t_comp(jj) = NaN; % Nan for open cells
            if t_filled(jj)
                t_comp(jj) = lego_inv(ismember(cell_list,t_coord(jj,:),'rows'),jj); % retrieve gface
            end
        end
        
        % intelligently choose/check shapes
        shp_ind = randperm(4); % randomize shape order
        rot_ang_ind = ang(randperm(4)); % randomize rotation
        
        flag1 = 0;
        while flag1 == 0
            c1 = 0;
            while c1 < 4 && flag1 == 0
                c1 = c1 + 1;
                shp = shp_ind(c1);
                c2 = 0;
                while c2 < 4 && flag1 == 0
                    c2 = c2 + 1;
                    rot_ang = rot_ang_ind(c2);
                    shape_cand = rot(shape(shp),shp,rot_ang); % candidate shape
                    if strcmp(rule,'strict')
                        tst = all(t_comp(~isnan(t_comp)) & shape_cand.lego(~isnan(t_comp))); % shape fits - conservative rule (few holes);
                    elseif strcmp(rule,'relaxed')
                        tst = all(t_comp(~isnan(t_comp)) == shape_cand.lego(~isnan(t_comp))); % shape fits - nice rule!
                    end
                    if tst
                        flag1 = 1;
                        ii = ii + 1;
                        cell_list(ii,:) = target;
                        shape_list(ii) = shape_cand;
                        lego(ii,:) = shape_cand.lego;
                        gface{ii} = find(shape_list(ii).lego);
%                       patch(shape_list(ii).coord(1,:) + target(1),shape_list(ii).coord(2,:) + target(2),'k')
                        for kk=gface{ii}
                            if kk == 1
                                target_buff = cat(1,target_buff,[target(1)-1 target(2)]);
                            elseif kk == 2
                                target_buff = cat(1,target_buff,[target(1) target(2)+1]);
                            elseif kk == 3
                                target_buff = cat(1,target_buff,[target(1)+1 target(2)]);
                            elseif kk == 4
                                target_buff = cat(1,target_buff,[target(1) target(2)-1]);
                            end
                        end
                    end
                end
            end
            flag1 = 1;
            target_buff(1,:) = []; % clear top of buffer
        end
    end
end

% check that the shape matches the desired size
if (max(cell_list(:,1)) - min(cell_list(:,1))) < (sz_r(1)-1) || (max(cell_list(:,2)) - min(cell_list(:,2))) < (sz_r(2)-1)
    gene = shape_proc(num,sz,rule,sym,draw,prnt)
    return
end

% add symmetry
if sym(1) && ~sym(2) % horizontal symmetry
    cell_list = cat(1,cell_list,[cell_list(:,1) (cell_list(:,2).*-1)-1]);
    for ii=1:size(shape_list,2)
        shape_list(ii+size(cell_list,1)/2) = flipp(shape_list(ii),180);
%         patch(shape_list(ii+size(cell_list,1)/2).coord(1,:) + cell_list(ii+size(cell_list,1)/2,1),shape_list(ii+size(cell_list,1)/2).coord(2,:)+ cell_list(ii+size(cell_list,1)/2,2),'k')
    end
elseif sym(2) && ~sym(1) % vertical symmetry
    cell_list = cat(1,cell_list,[(cell_list(:,1).*-1)-1 cell_list(:,2)]);
    for ii=1:size(shape_list,2)
        shape_list(ii+size(cell_list,1)/2) = flipp(shape_list(ii),90);
%         patch(shape_list(ii+size(cell_list,1)/2).coord(1,:) + cell_list(ii+size(cell_list,1)/2,1),shape_list(ii+size(cell_list,1)/2).coord(2,:)+ cell_list(ii+size(cell_list,1)/2,2),'k')
    end
elseif sym(1) && sym(2) % dual symmetry
    cell_list = cat(1,cell_list,[cell_list(:,1) (cell_list(:,2).*-1)-1]);
    for ii=1:size(shape_list,2)
        shape_list(ii+size(cell_list,1)/2) = flipp(shape_list(ii),180);
%         patch(shape_list(ii+size(cell_list,1)/2).coord(1,:) + cell_list(ii+size(cell_list,1)/2,1),shape_list(ii+size(cell_list,1)/2).coord(2,:)+ cell_list(ii+size(cell_list,1)/2,2),'k')
    end
    cell_list = cat(1,cell_list,[(cell_list(:,1).*-1)-1 cell_list(:,2)]);
    for ii=1:size(shape_list,2)
        shape_list(ii+size(cell_list,1)/2) = flipp(shape_list(ii),90);
%         patch(shape_list(ii+size(cell_list,1)/2).coord(1,:) + cell_list(ii+size(cell_list,1)/2,1),shape_list(ii+size(cell_list,1)/2).coord(2,:)+ cell_list(ii+size(cell_list,1)/2,2),'k')
    end
end

if draw
    figure;
    axis equal
    axis off
    for ii=1:size(shape_list,2)
        patch(shape_list(ii).coord(1,:) + cell_list(ii,1),shape_list(ii).coord(2,:)+ cell_list(ii,2),'k')
    end
    if prnt
        set(gcf,'paperpositionmode','auto');
        print(['./FARM/' num2str(num) '_' rule '_box_' num2str(sz(1)) 'x' num2str(sz(2)) '_' num2str(sym(1)) '_' num2str(sym(2))],'-dpdf');
    end
elseif prnt
    figure('visible','off');
    axis equal
    axis off
    for ii=1:size(shape_list,2)
        patch(shape_list(ii).coord(1,:) + cell_list(ii,1),shape_list(ii).coord(2,:)+ cell_list(ii,2),'k')
    end
    set(gcf,'paperpositionmode','auto');
    print(['./FARM/' num2str(num) '_' rule '_box_' num2str(sz(1)) 'x' num2str(sz(2)) '_' num2str(sym(1)) '_' num2str(sym(2))],'-dpdf');
end

% convert to gene (a 2d representation of each shape/orientation)

gene = zeros(sz(1),sz(2),3);
for ii=1:size(shape_list,2)
    gene(cell_list(ii,1) - min(cell_list(:,1))+1,cell_list(ii,2) - min(cell_list(:,2))+1,:) = round([shape_list(ii).id; shape_list(ii).ori])'; % rounding to kill the numerical error
end

gene(gene == -1) = 0; % converting into a binary-like representation

for ii=1:sz(1)
    for jj=1:sz(2)
        gene_bi(ii,jj) = bi2de([gene(ii,jj,2) gene(ii,jj,3)]); % convert binary digit to decimal
    end
end

gene = round(gene(:,:,1).*10 + gene_bi); % here we convert each shape (first digit) and orientation (second digit) to a unique identifier
gene = flipud(gene'); % re-orient to a cartesian representation
