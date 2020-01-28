function pixel_fun(num,sz)

% Genrate procedural (evolutionary?) shape generation

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

empt.coord = [0 0 0 0; 0 0 0 0];
empt.area = 0;
empt.lego = [1 1 1 1];

shape(1) = sqr;
shape(2) = tri;
shape(3) = cv;
shape(4) = cc;
% shape(5) = empt;

%% start procedure

f = figure('visible','off');
axis equal

target = [0 0]; % initial shape origin (bottom/lefthand)
ang = [0 90 180 270];

flag1 = 0;
flag2 = 0;

ii = 1;
shp = randsample(1:4,1);
rot_ang = randsample(ang,1);

cell_list(ii,:) = target;
shape_list(ii) = rot(shape(shp),shp,rot_ang);
gface{ii} = find(shape_list(ii).lego);
lego(ii,:) = shape_list(ii).lego;

patch(shape_list(ii).coord(1,:) + target(1),shape_list(ii).coord(2,:) + target(2),'k')

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

% Main drawing loop after initialization
while ~isempty(target_buff)
    target_buff = target_buff(randperm(size(target_buff,1)),:); % randomize order for even growth
    target_buff(ismember(target_buff,cell_list,'rows'),:) = []; % remove full possible targets
    if ~isempty(target_buff)
        target = target_buff(1,:); % choose next position
    else
        pixel_fun(num,sz)
        %         axis equal
        %         axis off
        %         set(gcf,'paperpositionmode','auto');
        %         print(['./FARM/' num2str(num)],'-dpng');
        %         close all
        return
    end
    % check that max size is not exceeded
    while (max([cell_list(:,1); target(1)]) - min([cell_list(:,1); target(1)])) >= sz(1) || (max([cell_list(:,2); target(2)]) - min([cell_list(:,2); target(2)])) >= sz(2)
        target_buff(1,:) = [];
        if isempty(target_buff)
            pixel_fun(num,sz)
            %             axis equal
            %             axis off
            %
            %             set(gcf,'paperpositionmode','auto');
            %             print(['./FARM/' num2str(num)],'-dpng');
            %             close all
            return
        end
        target = target_buff(1,:); % choose next position
    end
    
    % query compatibility conditions of new target cell (so conditions
    % outside)
    t_coord = [target(1) - 1 target(2) + 0; target(1) + 0 target(2) + 1; target(1) + 1 target(2) + 0; target(1) + 0 target(2) - 1];
    t_filled = ismember(t_coord,cell_list,'rows')'; % find filled cells (L/T/R/B)
    
    lego_inv = lego(:,[3 4 1 2]); % Needs to be re-ordered to R/B/L/T
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
                % TRICKY PART
                %               if all((t_comp(~isnan(t_comp)) - shape_cand.lego(~isnan(t_comp))) >= 0) % shape fits - more liberal rule - kinda cool!
                if all(t_comp(~isnan(t_comp)) & shape_cand.lego(~isnan(t_comp))) % shape fits
                    flag1 = 1;
                    flag2 = 1;
                    ii = ii + 1;
                    cell_list(ii,:) = target;
                    shape_list(ii) = shape_cand;
                    lego(ii,:) = shape_cand.lego;
                    gface{ii} = find(shape_list(ii).lego);
                    % Draw
                    patch(shape_list(ii).coord(1,:) + target(1),shape_list(ii).coord(2,:) + target(2),'k')
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

if (max([cell_list(:,1); target(1)]) - min([cell_list(:,1); target(1)])) >= sz(1) || (max([cell_list(:,2); target(2)]) - min([cell_list(:,2); target(2)])) >= sz(2)
    pixel_fun(num,sz)
    return
end

axis equal
axis off

set(gcf,'paperpositionmode','auto');
print(['./FARM/' num2str(num)],'-dpng');
close all

