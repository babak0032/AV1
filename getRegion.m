function [ ordered_mask, regions] = getRegion( Imwork, thresh, normalisation)
%getRegion Summary of this function goes here
%   Detailed explanation goes here
  Imback = double(imread('DATA1/bgframe.jpg','jpg'));

  % subtract background & select pixels with a big difference
  if normalisation
        rgb_sum_work = sum(Imwork, 3);
        sat_work = rgb_sum_work/3;  % and saturation
        rgb_sum_work(rgb_sum_work == 0) = 1;
        
        red_work = Imwork(:,:, 1)./rgb_sum_work;
        green_work = Imwork(:,:,2)./rgb_sum_work;
      
        rgb_sum_bg = sum(Imback, 3);
        sat_bg = rgb_sum_bg/3;  % and saturation
        rgb_sum_bg(rgb_sum_bg == 0) = 1;
        red_bg = Imwork(:,:,1)./rgb_sum_bg;
        green_bg = Imwork(:,:,2)./rgb_sum_bg;          
        
        
        fore = (abs(red_work-red_bg) > thresh) ...
        | (abs(green_work - green_bg) > thresh) ...
        | (abs(sat_work - sat_bg) > thresh);
  else
        fore = (abs(Imwork(:,:,1)-Imback(:,:,1)) > thresh) ...
        | (abs(Imwork(:,:,2) - Imback(:,:,2)) > thresh) ...
        | (abs(Imwork(:,:,3) - Imback(:,:,3)) > thresh);
  end

  % % Remove regions that are small
  foremm = bwareaopen(fore, 100);
  foremm = bwmorph(foremm,'close');
  foremm = bwmorph(foremm,'majority');
  foremm = bwareaopen(foremm, 50);
  
  
  foremm = bwlabel(foremm,4);
  mask = imfill(foremm,'holes');
  
  regions = regionprops(mask, 'Area', 'Centroid'); %Gets Area, centriod, and the bounding box
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
  ordered_mask = zeros(size(mask));
  for i = 1 : N
    ordered_mask(mask == id(i)) = i;
  end
  
  
end

