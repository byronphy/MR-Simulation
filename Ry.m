function Ry = Ry(angle) 
% rotation matrix about the y-axis

Ry = [cos(angle), 0, sin(angle);...
      0, 1, 0;...
      -sin(angle), 0, cos(angle)];

end