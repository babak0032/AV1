%Some varibales:
time_pause = [0.25,0.25];
Dancers = 4;
avg_colour = zeros(Dancers, 3, 256);
Ignore_pos = [121,222; 150,223];
next_dancer = 1;


% Start image
Imstart = imread(strcat('DATA1/frame110.jpg'),'jpg');

% show start image
% figure(1)
% imshow(Imstart)
% pause(1)

[mask, candidate_regions] = getRegion( double(Imstart), 10, 0);
imshow(mask)
pause(1)


for region = 1 : size(candidate_regions)
    if ~any(pdist2(Ignore_pos, [candidate_regions(region).Centroid(1), candidate_regions(region).Centroid(2)]) < 50)
        
        col = mask == region;
        imshow(col)
        pause(time_pause(1))
        
        redChannel = Imstart(:, :, 1);
        greenChannel = Imstart(:, :, 2);
        blueChannel = Imstart(:, :, 3);
        
        [redhist, ~] = imhist(redChannel(col), 256);
%         imhist(redChannel(col), 256)
%         pause(1)
        [greenhist, ~] = imhist(greenChannel(col), 256);
        [bluehist, ~] = imhist(blueChannel(col), 256);
        
        avg_colour(next_dancer,:,:) = [redhist, greenhist, bluehist]';
        
        next_dancer = next_dancer + 1;

    end
    if next_dancer > Dancers
        break
    end
end




% Overlapping image
Imwork = imread(strcat('DATA1/frame', int2str(19 + 109), '.jpg'),'jpg');

% Show the image
figure(2)
imshow(Imwork)
pause(1)

% Show the mask
[mask, candidate_regions] = getRegion( double(Imwork), 10, 0);
imshow(mask)
pause(1)


% split the region
[mask, props] = getRegion( double(Imwork), 10, 0);
cartProd = [2,2,1,1]


% First test check if region is occupied twice and alter mask
regions = unique(cartProd(1,:));

% If size of it is not 4 (some of them got aassigned twice)
if ~(size([regions],2) == Dancers)
    disp('Something')
    
    % size(cartProd(1,:)) % 1     4
    
    % Go through every tracker
    for track = 2 : size(cartProd(1,:),2)
        
        track
        
        % If a region a got picked twice
        if sum(cartProd(1,:) == cartProd(1,track)) > 1
            disp('Splitting')
            disp( cartProd(1,track) )
            nr = sum(cartProd(1,:) == cartProd(1,track)) % A bug here
            [mask, new_min, new_max] = split_region(mask, cartProd(1,track), 2);

        end
    end
    imshow(label2rgb(mask))
    pause(2)
end

disp('split over!')

% The colour test

% get the colour chennels
redChannel = Imwork(:, :, 1);
greenChannel = Imwork(:, :, 2);
blueChannel = Imwork(:, :, 3);
color_pass = 1;

size(candidate_regions);

% this matrix, for every dancer, store the colour difference of that dancer
% to every region
col_diffs = zeros(Dancers, size(candidate_regions,1));
min_ids = zeros(Dancers);

% for every dancer and every region
for d = 1 : Dancers
    for reg = 1 : size(candidate_regions,1)
        col = mask == reg;
        imshow(col)
        pause(time_pause(1))

        [redhist, ~] = imhist(redChannel(col), 256);
        [greenhist, ~] = imhist(greenChannel(col), 256);
        [bluehist, ~] = imhist(blueChannel(col), 256);
        
        % size(redhist) % (256, 1)

        % Store the difference of every bin and calculate the vector
        col_diffs(d,reg) = norm(redhist' - squeeze(avg_colour(d,1,:))') + ...
            norm(greenhist' - squeeze(avg_colour(d,2,:))') + ...
            norm(bluehist' - squeeze(avg_colour(d,3,:))');
    end
    [~,min_id] = min(col_diffs(d,:));
    disp(min_id)
    min_ids(d) = min_id;
end
col_diffs


pause(3)
for r = 1 : Dancers
    col = mask == min_ids(r);
    imshow(col)
    pause(1)

    [redhist, ~] = imhist(redChannel(col), 256);
    [greenhist, ~] = imhist(greenChannel(col), 256);
    [bluehist, ~] = imhist(blueChannel(col), 256);

    col_shift = norm(redhist' - squeeze(avg_colour(r,1,:))') + ...
            norm(greenhist' - squeeze(avg_colour(r,2,:))') + ...
            norm(bluehist' - squeeze(avg_colour(r,3,:))');
    col_shift;
    color_pass = color_pass & (col_shift < 160);
end

if ~color_pass
    disp('COLORS');
    continue
end


% For later


% Overlapping image
Imwork = imread(strcat('DATA1/frame', int2str(19 + 109), '.jpg'),'jpg');

% Show the image
figure(2)
imshow(Imwork)
pause(1)

% Show the mask
[mask, properties] = getRegion( double(Imwork), 10, 0);
imshow(mask)
pause(1)