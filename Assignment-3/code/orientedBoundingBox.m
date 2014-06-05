function [bb, elong] = orientedBoundingBox(pos, dir)
% Angle of group direction
v = -atan2(dir(2),dir(1));
% Rotation matrix
R = [cos(v), -sin(v); sin(v), cos(v)];
rotPos = zeros(size(pos));
% Rotate positions
for i=1:size(pos,2)
	rotPos(:,i) = R*pos(:,i);
end
[minX,minXidx] = min(rotPos(1,:));
[maxX,maxXidx] = max(rotPos(1,:));
[minY,minYidx] = min(rotPos(2,:));
[maxY,maxYidx] = max(rotPos(2,:));

%minX = pos(1,minXidx);
%maxX = pos(1,maxXidx);
%minY = pos(2,minYidx);
%maxY = pos(2,maxYidx);
v = -v;
R = [cos(v), -sin(v); sin(v), cos(v)];
bb(:,1) = R*[minX; minY];
bb(:,2) = R*[minX; maxY];
bb(:,3) = R*[maxX; maxY];
bb(:,4) = R*[maxX; minY];
elong = (maxX-minX)/(maxY-minY);

end