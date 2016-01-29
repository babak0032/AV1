
% compute the background image

Imback = double(imread('DATA1/bgframe.jpg','jpg'));
[MR,MC,Dim] = size(Imback);

% loop over all images
fig1=1;
fig2=1;
fig15=1;
fig3=1;
fig4=1;
for i = 110 : 110
  % load image
  Im = (imread(['DATA1/frame', int2str(i), '.jpg'],'jpg')); 

  if fig1 > 0
    figure(fig1)
    clf
    imshow(Im)
    pause(1)
  end
  Imwork = double(Im);

  %extract ball
  [cc(i),cr(i),radius,flag]=extractball(Imwork,Imback,fig1,fig2,fig3,fig15,i);
  if flag==0
    continue
  end

  if fig1 > 0
    figure(fig1)
    hold on
    for c = -0.97*radius: radius/20 : 0.97*radius
      r = sqrt(radius^2-c^2);
      plot(cc(i)+c,cr(i)+r,'g.')
      plot(cc(i)+c,cr(i)-r,'g.')
    end
    %eval(['saveas(gcf,''TRACK/trk',int2str(i-1),'.jpg'',''jpg'')']);  
  end

      pause(1.3)
end

% show positions
if fig4 > 0
  figure(fig4)
  hold on
  clf
  plot(cc,'r*')
  plot(cr,'g*')
end
