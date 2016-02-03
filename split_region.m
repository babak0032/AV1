function [ mask, indices ] = split_region( mask, region, dancers)
% split_region( mask, region, dancers)
%           mask - image containing the regions
%           region - index of the region to be split
%           dancers [Dancers, 2] - Dancers, their estimated coordinates
     new_index = max(max(mask)) - 1;
     dancers(2:end,:);
     [r, c] = find(mask == region);
     for pixel = 1 : size(r,1)
         new_i = knnsearch(dancers(:,:), [r(pixel), c(pixel)]);
         mask(r(pixel), c(pixel)) = new_i + new_index; 
     end
     indices = [region, new_index + 2 : new_index + size(dancers,1)];
end

