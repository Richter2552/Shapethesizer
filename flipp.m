function [output] = flipp(input,theta)

% theta = 180 for horizontal and 90 for vertical

theta_rad = deg2rad(theta);

output = input;
output.coord = [cos(2*theta_rad) sin(2*theta_rad); sin(2*theta_rad) -cos(2*theta_rad)] * input.coord;
output.coord = output.coord - min(output.coord,[],2);

if theta == 180
    output.lego = permute(input.lego,[1 4 3 2]);
elseif theta == 90
    output.lego = permute(input.lego,[3 2 1 4]);
end


if output.id > 1
    output.ori = [cos(2*theta_rad) sin(2*theta_rad); sin(2*theta_rad) -cos(2*theta_rad)] * output.ori; % reflection of orientation vector
else
    output.ori = [1; 1];
end
