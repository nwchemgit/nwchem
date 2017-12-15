function [target, time_stamps] = md_compress(trace,sampling_rate)

    % reformat data
    [x, y, z] = size(trace);
    t2 = reshape(trace, [x, y*z]);
    %% SVD Analysis
    [c,s,~] = svds(t2,x);


    %% project to low dimensions
    n_dim = 2;
    ps = c(:,1:n_dim) * s(1:n_dim, 1:n_dim);
    psd = ps(1:(x-1),1:n_dim) - ps(2:x, 1:n_dim);
    psdu = sum(abs(psd'))
    prob_dist = (psdu / sum(psdu));


    %% Sub-sampling
    total = 10^-10 + prob_dist(1);
    sampling_entries = int32(x * sampling_rate);
    target = zeros(sampling_entries, y, z);
    time_stamps = zeros(sampling_entries,1);
    %output initialization: the first one should be added always
    target(1,:,:) = trace(1,:,:);
    time_stamps(1) = 1;
    t_idx = 2;
    for i=2:x
        total = total + prob_dist(i-1);
        if(total > 1 / double(sampling_entries - 1))
            target(t_idx,:,:) = trace(i,:,:);
            time_stamps(t_idx) = i;
            total = total - 1 / double(sampling_entries);
            t_idx = t_idx + 1;
        end
    end

end