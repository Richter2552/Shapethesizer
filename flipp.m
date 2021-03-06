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
    output.ori = round([cos(2*theta_rad) sin(2*theta_rad); sin(2*theta_rad) -cos(2*theta_rad)] * input.ori); % reflection of orientation vector
    % conversion of vector to a handy notation (0,1,2,3) = 0,90,180,270 degrees
    ori_id = mod(atan2(output.ori(2),output.ori(1)),2*pi);
    output.ori_id = ori_id/(pi/2) - 0.5;
else
    output.ori = input.ori;
    output.ori_id = 0;
end
