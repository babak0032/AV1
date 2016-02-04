function [ mask, new_index ] = split_region( mask, region, props, n)
% split_region( mask, region, dancers)
%           mask - image containing the regions
%           region - index of the region to be split
%           dancers [Dancers, 2] - Dancers, their estimated coordinates
%      new_index = max(max(mask));
%      figure(2)
%      imshow(mask==region)
%      pause(1)
%      
%      [r, c] = find(mask == region);
%      for pixel = 1 : size(r,1)
%          mask(r(pixel), c(pixel)) = knnsearch(dancers(:,:), [c(pixel), r(pixel)]) + new_index;
%      end
%      imshow(label2rgb(mask))     
%      hold on
%      for dancer = 1 : 2
%          plot(dancers(dancer,1), dancers(dancer,2), 'r.')
%      end
%      pause(5)
%      hold off
    new_index = max(max(mask));
    
    dir = [cos(props(region).Orientation), sin(props(region).Orientation)];
    
    midpoints = zeros(n - 1, 2);
    midpoints(1,:) = [props(region).Centroid(1) - cos(props(region).Orientation) * props(region).MajorAxisLength * (n-2)/2/n, props(region).Centroid(2) - sin(props(region).Orientation) * props(region).MajorAxisLength * (n-2)/2/n];
    for i = 2 : n - 1
        midpoints(i,:) = [midpoints(i-1,1) + cos(props(region).Orientation) * props(region).MajorAxisLength /n, midpoints(i-1,2) + sin(props(region).Orientation) * props(region).MajorAxisLength /n];
    end
    
    
    [r, c] = find(mask == region);
    mask(mask == region) = new_index + n;
    for pixel = 1 : size(r,1)
        for i = 1 : size(midpoints,1)
            if ([midpoints(i,2), midpoints(i,1)] - [r(pixel), c(pixel)]) * dir' < 0
                mask( r(pixel), c(pixel)) = new_index + i;
                break
            end
        end
    end
end

