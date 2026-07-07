function [sort_edge,sort_point] = splitBoundaryEdges(B,point)
% SPLITBOUNDARYEDGES Split a boundary loop into four ordered edge segments.
% Author: Joey.Chow
location = [];
for k = 1:length(point)
    for i = 1:length(B)
        if B(i,1) == point(k)
            location = [location;i];
        end
    end
end
[sort_location,index] = sort(location);
sort_point = point(index);
for i = 1:length(sort_location)-1
    edge{i} = B(sort_location(i,1):sort_location(i+1,1),1);
end
edge{length(sort_location)} = [B(sort_location(end,1):end,1);B(1:sort_location(1,1),1)];
for i = 1:length(point)-1
    for j = 1:length(edge)
        temp = edge{j};
        if sum(ismember([point(i);point(i+1)],[temp(1);temp(end)])) == 2
            sort_edge{i} = temp;
        end
    end
end
for j = 1:length(edge)
    temp = edge{j};
    if sum(ismember([point(end);point(1)],[temp(1);temp(end)])) == 2
        sort_edge{length(point)} = temp;
    end
end
end
