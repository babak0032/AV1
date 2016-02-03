function [ regions ] = regions_sorted( mask )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
  regions = regionprops(mask, ['basic']); %Gets Area, centriod, and the bounding box
  [N,~] = size(regions);
  if N < 1
    return   
  end

  % do bubble sort (large to small) on regions in case there are more than 1
  id = zeros(N);
  for i = 1 : N
    id(i) = i;
  end
  for i = 1 : N-1
    for j = i+1 : N
      if regions(i).Area < regions(j).Area
        tmp = regions(i);
        regions(i) = regions(j);
        regions(j) = tmp;
        tmp = id(i);
        id(i) = id(j);
        id(j) = tmp;
      end
    end
  end

end

