function dummy_answer = check_track(num, path)
Colours = cellstr(['g.'; 'b.'; 'c.'; 'r.'; 'g.'; 'b.']);

Im = imread(['DATA1/frame', int2str(num), '.jpg'],'jpg');
track_num = 50;
imshow(Im)
pause(1)
hold on

index = 1;

% for i = num:num+10
%   for region = 1:2
%     char(Colours(index));
%     radius = 10;
%     plot(path(region,i-109,1),path(region,i-109,2),char(Colours(index)))
%     for c = -0.99*radius: radius/10 : 0.99*radius
%       r = sqrt(radius^2-c^2);
%       plot(path(region,i-109,1) + c, path(region,i-109,2) + r, char(Colours(index)))
%       plot(path(region,i-109,1) + c, path(region,i-109,2) - r, char(Colours(index)))
%     end
%     index = index + 1;
% 
%   end
%   index = 1;
% end

plot(squeeze(path(1,num-109,1)), squeeze(path(1,num-109,2)), 'r*')
plot(squeeze(path(1,num-109:num+track_num-109,1)), squeeze(path(1,num-109:num+track_num-109,2)), 'r')
plot(squeeze(path(2,num-109,1)), squeeze(path(2,num-109,2)), 'b*')
plot(squeeze(path(2,num-109:num+track_num-109,1)), squeeze(path(2,num-109:num+track_num-109,2)), 'b')
% plot(squeeze(path(3,:,1)), squeeze(path(3,:,2)), 'g')
% plot(squeeze(path(4,:,1)), squeeze(path(4,:,2)), 'c')
pause(5)

end