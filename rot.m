function [output] = rot(input,ind,theta)

% note this can only shift by intervals of 90 degrees

theta_rad = deg2rad(theta);

output = input;
output.coord = [cos(theta_rad) -sin(theta_rad); sin(theta_rad) cos(theta_rad)]*input.coord;
output.coord = output.coord - min(output.coord,[],2);

output.lego = circshift(input.lego,-round(theta_rad/(pi/2)));

output.id = ind;

if output.id > 1
    output.ori = [cos(theta_rad) -sin(theta_rad); sin(theta_rad) cos(theta_rad)]*[1; 1]; % rotation of origin vector
else
    output.ori = [1; 1];
end

