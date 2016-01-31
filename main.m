
Imback = double(imread('DATA1/bgframe.jpg','jpg'));
[MR,MC,Dim] = size(Imback)

fig1=0;
fig2=0;
fig15=0;
fig3=0;
fig4=1;

NCON = 100;

% Init state vectors and weights and predicted estimates
x = zeros(NCON,4);
P = zeros(NCON,4,4);
weights = zeros(NCON,1);

for i = 1:NCON
  OP(i,1,1) = 100;
  OP(i,2,2) = 100;
  OP(i,3,3) = 100;
  OP(i,4,4) = 100;
end

for k = 1 : NCON
  x(k,:) = [floor(MC*rand(1)),floor(MR*rand(1)),0,0]';
  weights(k)=1/NCON;
end


for i=110:169
    Im = (imread(['DATA1/frame',int2str(i), '.jpg'],'jpg'));
    if fig1 > 0
      figure(fig1)
      clf
      imshow(Im)
    end
    Imwork = double(Im);
    
    [cc,cr,radius,flag]=extractball(Imwork,Imback,fig1,fig2,fig3,fig15,i);
    z = [cc,cr,radius]; % Observation
    tracks1(i-109) = cc;
    tracks2(i-109) = cr;
    
    if fig1 > 0
      figure(fig1)
      hold on
      for c = -0.99*radius: radius/10 : 0.99*radius
        r = sqrt(radius^2-c^2);
        plot(cc+c,cr+r,'g.')
        plot(cc+c,cr-r,'g.')
      end
    end
    
    [x,P,weights] = condense_function(x,P,weights,z);  
    
    top = find(weights == max(weights));
    weights(top);
    x(top,1);
    x;
    
    % tracks1(i-109) = x(top,1);
    % tracks2(i-109) = x(top,2);
    
    if fig1 > 0
      figure(fig1)
      hold on
      for c = -0.99*radius: radius/10 : 0.99*radius
        r = sqrt(radius^2-c^2);
    %    plot(x(top,i,1)+c,x(top,i,2)+r+1,'b.')
    %    plot(x(top,i,1)+c,x(top,i,2)+r,'y.')
        plot(x(top,1)+c,x(top,2)+r,'r.')
        plot(x(top,1)+c,x(top,2)-r,'r.')
    %    plot(x(top,i,1)+c,x(top,i,2)-r,'y.')
    %    plot(x(top,i,1)+c,x(top,i,2)-r-1,'b.')
      end
    %  eval(['saveas(gcf,''COND/cond',int2str(i-1),'.jpg'',''jpg'')']);  
    end
    
end

tracks1;

if fig4 > 0
  figure(fig4)
  hold on
  clf
  plot(tracks1, tracks2, 'r*')
end