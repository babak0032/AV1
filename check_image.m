function dummy_answer = check_image(num, path, actual, split_mode)

Colours = cellstr(['g.'; 'b.'; 'c.'; 'r.'; 'g.'; 'b.']);
Colours2 = cellstr(['g*'; 'b*'; 'c*'; 'r*'; 'r*']);
Ignore_pos = [121,222; 150,223];
index = 1;


Imstart = imread(['DATA1/frame', int2str(num), '.jpg'],'jpg');
figure(1)
imshow(Imstart)
pause(1)
hold on

figure(2)

[mask_test, candidate_regions_test] = getRegion( double(Imstart), 12, 0);

imshow(mask_test)
pause(1)

if split_mode
    [mask_test, new_min, new_max] = split_region(mask_test, 1, 2);
    [mask_test, new_min, new_max] = split_region(mask_test, 2, 2);
    figure(3)
    imshow(label2rgb(mask_test))
    pause(2)
    
    candidate_regions_test = regionprops(mask_test, 'Area', 'Centroid'); %Gets Area, centriod, and the bounding box
end



figure(4)
imshow(Imstart)
pause(1)
hold on

for region = 1:size(candidate_regions_test,1)
    if ~any(pdist2(Ignore_pos, [candidate_regions_test(region).Centroid(1), candidate_regions_test(region).Centroid(2)]) < 170)
        radius = sqrt(candidate_regions_test(region).Area / pi);
        plot(candidate_regions_test(region).Centroid(1),candidate_regions_test(region).Centroid(2),char(Colours(index)))
        for c = -0.99*radius: radius/10 : 0.99*radius
          r = sqrt(radius^2-c^2);
          plot(candidate_regions_test(region).Centroid(1) + c, candidate_regions_test(region).Centroid(2) + r, char(Colours(index)))
          plot(candidate_regions_test(region).Centroid(1) + c, candidate_regions_test(region).Centroid(2) - r, char(Colours(index)))
        end
        index = index + 1;
    end
end
pause(2)

% for region = 1:4
% 
%     radius = 10;
%     plot(path(region,num-109,1),path(region,num-109,2),char(Colours(index)))
%     for c = -0.99*radius: radius/10 : 0.99*radius
%       r = sqrt(radius^2-c^2);
%       plot(path(region,num-109,1) + c, path(region,num-109,2) + r, char(Colours(index)))
%       plot(path(region,num-109,1) + c, path(region,num-109,2) - r, char(Colours(index)))
%     end
%     index = index + 1;
% 
% end
% pause(2)

% hold on
% 

% Draw your observations
% for i=1:4
%    plot(path(i,num-109,1),path(i,num-109,2),'y*')
% end
% pause(2)

% Draw actual positions
% for i=1:4
%    plot(path(i,num-109,1),actual(i,num-109,2),char(Colours2(i)))
% end
% pause(2)

% figure(2)
% imshow(mask_test)
% pause(0.1)

end