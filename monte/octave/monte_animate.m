function monte_animate(t,X,Z,u,uscale,fps)

dt = t(2)-t(1);

xorder = 0;
xorderold = 0;
vorder = 0;
vorderold = 0;
sorder = 0;
sorderold = 0;

variances = var(X);
skews = skewness(X);

for i=2:length(t)
     subplot(3,1,1)
     
     scatter(X(:,i),Z(:,i))
     % Plot the mean velocity line in red.
     
     % xmean = 0 as long as we use zero-mean flow.
     xmean = 0;
     sample = linspace(0,1,101);
     hold on
%     plot([xmean,xmean],[0,1],'r')
     plot(u(sample)*t(i),sample,'r')
     hold off
     
     % Be clever about the axes.
     xorder = ceil(log10(max(abs(X(:,i)))) / log10(4));
     xorder = max(xorderold,xorder);
     xorderold = xorder;
     
     bound = 4^xorder;
     axis([-bound, bound, 0, 1])
     
     % Make an animated title.
     mytitle = strcat('t = ',num2str(t(i)));
     title(mytitle)
     xlabel('x'), ylabel('z')
     
     % Variance
     subplot(3,1,2)
     
     vorder = ceil(log10(max(abs(variances(:,i)))) / log10(4));
     vorder = max(vorderold,vorder);
     vorderold = vorder;
     
     plot(t(1:i), variances(1:i))
     axis([0, 4^ceil( log10(t(i))/log10(4) ), 0, 4^vorder])
     xlabel('t'), ylabel('Variance')
     grid on
     
     % Skewness
     subplot(3,1,3)
     
     sorder = ceil(log10(max(abs(skews(:,i)))) / log10(4));
     sorder = max(sorderold,sorder);
     sorderold = sorder;
     
     plot(t(1:i), skews(1:i))
     axis([0, 4^ceil( log10(t(i))/log10(4) ), 0, max(skews(2:i))])
     xlabel('t'), ylabel('Variance')
     grid on
     
     pause(1/fps)
end

end
