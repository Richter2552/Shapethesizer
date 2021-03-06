function [output] = rot(input,ind,theta)

% note this can only shift by intervals of 90 degrees

theta_rad = deg2rad(theta);

output = input;
output.coord = [cos(theta_rad) -sin(theta_rad); sin(theta_rad) cos(theta_rad)]*input.coord;
output.coord = output.coord - min(output.coord,[],2);

output.lego = circshift(input.lego,-round(theta_rad/(pi/2)));

output.id = ind;

if output.id > 1
    output.ori = round([cos(theta_rad) -sin(theta_rad); sin(theta_rad) cos(theta_rad)]*[1; 1]); % rotation of origin vector
    % conversion of vector to a handy notation (0,1,2,3) = 0,90,180,270 degrees
    ori_id = mod(atan2(output.ori(2),output.ori(1)),2*pi);
    output.ori_id = ori_id/(pi/2) - 0.5;
else
    output.ori = [1; 1];
    output.ori_id = 0;
end


