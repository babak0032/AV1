function [ mask, new_min, new_max ] = split_region( mask, region, n)
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

    props = regionprops(mask == region, 'Centroid','Orientation', 'MajorAxisLength');
    dir = [cos(degtorad(props.Orientation)); sin(-degtorad(props.Orientation))];
    
    midpoints = zeros(n - 1, 2);
    midpoints(1,:) = [props.Centroid(1) - cos(degtorad(props.Orientation)) * props.MajorAxisLength * (n-2)/2/n, props.Centroid(2) - sin(-degtorad(props.Orientation)) * props.MajorAxisLength * (n-2)/2/n];
    for i = 2 : n - 1
        midpoints(i,:) = [midpoints(i-1,1) + cos(degtorad(props.Orientation)) * props.MajorAxisLength / n, midpoints(i-1,2) + sin(-degtorad(props.Orientation)) * props.MajorAxisLength / n];
    end
    
    
    [r, c] = find(mask == region);
    mask(mask == region) = new_index + n;
    for pixel = 1 : size(r,1)
        for i = 1 : size(midpoints,1)
            if ([ c(pixel) - midpoints(i,1), r(pixel) - midpoints(i,2)] * dir) < 0
                mask( r(pixel), c(pixel)) = new_index + i;
                
                
                break
            end
        end
    end
    new_min = new_index + 1;
    new_max = new_index + n;

end

