function [] = plotseries(ffrst,flast,fskip)


% Output:  g0000.jpg, g0001.jpg, g0002.jpg, ...
% Can join as mpeg, e.g.
%    avconv -r 25 -i g%04d.jpg -b 1800k out.mp4

ct = 0;

for fl = ffrst:fskip:flast

   % RESET, SHOW/HIDE FIGURE
   close all;
%   figure('Position',[0 0 880 880]);
   figure('Position',[0 0 880 880],'Visible','Off');

   % PLOTTING 
   flname = sprintf('mat_vec%04i.cdf',fl);
   isosurf(flname,0.8);
   hold on;
   % ... add more to figure here


   % SAVE CURRENT FIGURE TO JPEG: g????.jpg
   flname = sprintf('g%04i.jpg', ct);
   print('-djpeg','-r150',flname);
   ct = ct + 1;

end 

end
