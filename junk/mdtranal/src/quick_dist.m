%% Set-up
% step size for timestamp
step = 100;
load data/nt-ice_md-normal_new.mat
% sampling rate 
sampling_rate = 0.1; % 10%
% reformat data
[x, y, z] = size(trace);
t2 = reshape(trace, [x, y*z]);
%% SVD Analysis
[c,s,l] = svd(t2);

% Visualize the results in 2D
figure, plot(c(:,1), c(:,2))
hold on
for i =1:step:x
    text(c(i,1),c(i,2), int2str(i))
end

% Visualize the results in 3D
figure, plot3(c(:,1), c(:,2), c(:,3))
hold on
for i =1:step:x
    text(c(i,1),c(i,2),c(i,3), int2str(i))
end

%% project to low dimensions
n_dim = 2;
ps = c(:,1:n_dim) * s(1:n_dim, 1:n_dim);
%% visualize changes
psd = ps(1:(x-1),1:n_dim) - ps(2:x, 1:n_dim);
psdu = sum(abs(psd'));
prob_dist = (psdu / sum(psdu));
figure, plot(1:(x-1), sum(abs(psd),2));


%% Sub-sampling
total = 10^-10;
sampling_entries = int32(x * sampling_rate);
target = zeros(sampling_entries, y, z);
t_idx = 1;
for i=1:(x-1)
    total = total + prob_dist(i);
    if(total > 1 / double(sampling_entries))
        target(t_idx,:,:) = trace(i,:,:);
        total = total - 1 / double(sampling_entries);
        t_idx = t_idx + 1;
    end
end
% TODO: post processing of trimming may be required as t_idx may be ended
% earlier due to one entry's prob_dist(i) can be 2 or more times than 1 /
% sampling_entries


