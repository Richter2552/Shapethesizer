function gene_render(gene)

% gene_render contructs patch images from gene for export or further
% processing

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

sz = size(gene);
gene = flipud(gene); % re-orient the gene

shp_id = zeros(sz(1),sz(2));
ori_id = zeros(sz(1),sz(2));
for ii=1:sz(1)
    for jj=1:sz(2)
        shp_id(ii,jj) = floor(gene(ii,jj)/10); % get shape id
        ori_id(ii,jj) = mod(gene(ii,jj),10); % ori_id
    end
end

shp_id(shp_id == -1) = 0; % convert negative space to 0

h = figure

for ii=1:sz(1)
    for jj=1:sz(2)
        if shp_id(ii,jj) > 0
            theta_rad = ori_id(ii,jj) .* pi/2;
            shp_rot = [cos(theta_rad) -sin(theta_rad); sin(theta_rad) cos(theta_rad)] * shape(shp_id(ii,jj)).coord; % rotate coordinates
            shp_rot = shp_rot - min(shp_rot,[],2);
            patch(shp_rot(1,:) + jj, shp_rot(2,:) + ii,'k')
        end
    end
end

end


