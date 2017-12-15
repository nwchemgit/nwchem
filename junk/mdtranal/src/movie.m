%% Set-up
% step size for timestamp
step = 100;
load data/nt-ice_md-normal_new.mat
% reformat data
[x, y, z] = size(trace);
t3 = reshape(trace, [x*y,z]);
% svd
[u, s, v] = svds(t3,2);
t4 = reshape(u,[x,y,2]);

% set 
xmin = min(min(t4(:,:,1)));
xmax = max(max(t4(:,:,1)));

ymin = min(min(t4(:,:,2)));
ymax = max(max(t4(:,:,2)));


%% Generate movie clip.
writerObj = VideoWriter('out.avi');
writerObj.FrameRate = 30;
open(writerObj); 
fid = figure; 
for i=1:x      
    figure(fid); % Makes sure you use your desired frame.
    plot(t4(i,:,1),t4(i,:,2),'.');
    title(sprintf('Frame: %04d',i));
    % ensure ranges
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end

close(writerObj); % Saves the movie.
