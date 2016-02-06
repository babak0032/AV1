function answer = evaluation(regions, track_pos)
% Plot every thing
% regions: all 4 regions (no order of dancers) (frames, regions, state_vector)
% track_pos:(dancers, particles, frames, state_vector) : (4,200,210,6)

% vector of number of correct trackings for each dancer
matchings = zeros(4);
dist_stat = zeros(4,2);

load('DATA1/positions1.mat')
answer = 1;

for frame=1:210
  % plot the ground truth positions (4,210,2)
  % plot(squeeze(positions(dancer,:,1)), squeeze(positions(dancer,:,2)), 'r*')
  % pause(2)
  
  region_pos = squeeze(regions(frame,:,:));
  size(region_pos);
  
  region_pos_x = region_pos(:,1);
  region_pos_y = region_pos(:,2);
  size([region_pos_x, region_pos_y]);
  region_pos = [region_pos_x, region_pos_y];
  
  
  distances = pdist2(squeeze(positions(:,frame,:)), region_pos, 'euclidean');
  size(distances);
  
  matchings = sum((distances < 10),2)>0;
  matchings
  
  %dist_stat(dancer,:) = [mean(distances), std(distances)]'
  
  
  
end



end