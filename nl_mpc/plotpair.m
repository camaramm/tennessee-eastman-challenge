function plotpair(t,y1,y2)

% Plots pairs from y1 and y2
%
%  plotpair(t,y1,y2)

if min(size(t)) ~= 1
   error('T must be a vector')
end
M=length(t);

if any(size(y1)-size(y2))
   error('Y1 and Y2 must be same size')
end

[m1,n1]=size(y1);
if m1 ~= M
   error('# rows in Y1 must equal length of T')
end

if n1 == 1
   clg
   plot(t,[y1 y2])
   title(['Variable # 1'])
elseif n1 == 2
   clg
   for j=1:2
      subplot(210+j)
      plot(t,[y1(:,j) y2(:,j)])
      title(['Variable # ',int2str(j)])
   end
else
   for i=1:4:n1
      clg
      for j=0:min([3,n1-i])
         subplot(221+j)
         plot(t,[y1(:,i+j) y2(:,i+j)])
         title(['Variable # ',int2str(i+j)])
      end
      if i+j < n1, pause, end
   end
end
