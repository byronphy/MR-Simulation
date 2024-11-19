function Rx = Rx(angle) 
% rotation matrix about the x-axis

Rx = [1, 0, 0;...
      0, cos(angle), -sin(angle);...
      0, sin(angle), cos(angle)];

end