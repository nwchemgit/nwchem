%% Set-up
% step size for timestamp
step = 40;
load data//aimd_nt_dat_noc.mat
sampling_rates = [0.01 0.05 0.06 0.07 0.08];
%sampling_rates = [0.05];
[x, y, z] = size(trace);
t2 = reshape(trace, [x, y*z]);
%% SVD Analysis
[c,s,l] = svd(t2);

rates_ssg = zeros(size(sampling_rates));
rates_deg = zeros(size(sampling_rates)); % differential entropy of Gaussian
tss = cell(size(sampling_rates));
for i=1:size(sampling_rates,2)
    [target, tss{i}] = md_compress_acc(trace, sampling_rates(i));
    usamples = trace(int32(linspace(1,x,size(target,1))),:,:);
    sampling_rates(i)
    size(target,1)
    size(usamples, 1)
    
    assert(size(target,1) == size(usamples,1));
    cssg = eval_singular_gap(target);
    ussg = eval_singular_gap(usamples);
    rates_ssg(i) = cssg / ussg;
    cdeg = eval_diff_ent_gau(target);
    udeg = eval_diff_ent_gau(usamples);
    rates_deg(i) = cdeg / udeg;
    sprintf("sr:%f ssgr=%f/%f=%f degr=%f/%f=%f",sampling_rates(i), cssg, ussg, rates_ssg(i), cdeg, udeg, rates_deg(i))
    % Visualize the results in 2D
    f = figure, plot(c(:,1), c(:,2))
    hold on
    for j =1:step:x
        text(c(j,1),c(j,2), int2str(j))
    end
    
    plot(c(tss{i},1), c(tss{i},2), 'rx')
    plot(c(int32(linspace(1,x,size(target,1))),1), c(int32(linspace(1,x,size(target,1))), 2), 'bo')
    hold off
end 

%% Draw result plot
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create loglog
loglog(sampling_rates,rates_ssg,'Marker','o','LineStyle','-','Color',[1 0 0]);
loglog(sampling_rates,rates_deg,'Marker','*','LineStyle','-','Color',[0 0 1]);
% Create xlabel
xlabel({'Sampling Rate (log)'});

% Create ylabel
ylabel({'SSGR/DEGD'});

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',12,'XMinorTick','on','XScale','log','YMinorTick','on',...
    'YScale','linear');


