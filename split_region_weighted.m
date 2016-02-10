function [ mask, new_min, new_max ] = split_region_weighted( mask, region, weights)
    new_index = max(max(mask));
    n = size(weights, 2);
    weights = cumsum(weights);
    props = regionprops(mask == region, 'Centroid','Orientation', 'MajorAxisLength');
    dir = [cos(degtorad(props.Orientation)); sin(-degtorad(props.Orientation))];
    
    midpoints = zeros(n - 1, 2);
    end_point = [props.Centroid(1) - cos(degtorad(props.Orientation)) * props.MajorAxisLength / 2, props.Centroid(2) - sin(-degtorad(props.Orientation)) * props.MajorAxisLength /2];
    for i = 1 : n - 1
        midpoints(i,:) = [end_point(1) + cos(degtorad(props.Orientation)) * props.MajorAxisLength * weights(i), end_point(2) + sin(-degtorad(props.Orientation)) * props.MajorAxisLength * weights(i)];
        plot(midpoints(i,1), midpoints(i,2), 'r.')
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