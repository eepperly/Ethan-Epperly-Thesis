thesis_startup

s_list = 20:20:300;
num_trials = 1000;

hutchpp_errs = zeros(length(s_list), num_trials, 2);
hutchpp_iso_errs = zeros(length(s_list), num_trials, 2);
xtrace_errs = zeros(length(s_list), num_trials, 2);
xtrace_iso_errs = zeros(length(s_list), num_trials, 2);
xnystrace_errs = zeros(length(s_list), num_trials, 2);
xnystrace_iso_errs = zeros(length(s_list), num_trials, 2);

for i = 1:2
    if i == 1
        A = rand_with_evals(linspace(1,3,n));
    else
        A = rand_with_evals([ones(1,ceil(n/20)) 1e-3*ones(1,n-ceil(n/20))]);
    end
    trA = trace(A);
    
    for s_idx = 1:length(s_list)
        s = s_list(s_idx);
        for trial = 1:num_trials
            [s trial]
            hutchpp_errs(s_idx,trial,i) = abs(hutchpp(@(X) A*X,n,s)-trA) / trA;
            hutchpp_iso_errs(s_idx,trial,i) = abs(hutchpp_resphere(@(X) A*X,n,s)-trA) / trA;
            xtrace_errs(s_idx,trial,i) = abs(xtrace(@(X) A*X,n,s)-trA) / trA;
            xtrace_iso_errs(s_idx,trial,i) = abs(xtrace_resphere(@(X) A*X,n,s)-trA) / trA;
            xnystrace_errs(s_idx,trial,i) = abs(xnystrace(@(X) A*X,n,s)-trA) / trA;
            xnystrace_iso_errs(s_idx,trial,i) = abs(xnystrace_resphere(@(X) A*X,n,s)-trA) / trA;
        end
    end
end

%% Plot 1

figure(1)

quantiles = quantile(hutchpp_errs(:,:,1),0.5,2);
plot(s_list,quantiles,"Color",purple,"LineStyle","-","Marker","*","MarkerSize",10); hold on
quantiles = quantile(hutchpp_iso_errs(:,:,1),0.5,2);
plot(s_list,quantiles,"Color",darkpurple,"LineStyle","--","Marker","*","MarkerSize",10)

quantiles = quantile(xtrace_errs(:,:,1),0.5,2);
plot(s_list,quantiles,"Color",blue,"LineStyle","-","Marker","x","MarkerSize",10)
quantiles = quantile(xtrace_iso_errs(:,:,1),0.5,2);
plot(s_list,quantiles,"Color",darkblue,"LineStyle","--","Marker","x","MarkerSize",10)

quantiles = quantile(xnystrace_errs(:,:,1),0.5,2);
plot(s_list,quantiles,"Color",orange,"LineStyle","-","Marker","o","MarkerSize",10)
quantiles = quantile(xnystrace_iso_errs(:,:,1),0.5,2);
plot(s_list,quantiles,"Color",darkorange,"LineStyle","--","Marker","o","MarkerSize",10)

set(gca,"YScale","log")

xlabel("Number of matvecs $s$")
ylabel("Relative error")

exportgraphics(gcf,"../figs/fig_14_2_a.png")
saveas(gcf,"../figs/fig_14_2_a.fig")

%% Plot 2

figure(2)

quantiles = quantile(hutchpp_errs(:,:,2),[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"LineStyle","-","Marker","*","MarkerSize",10)
quantiles = quantile(hutchpp_iso_errs(:,:,2),[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),darkpurple,"LineStyle","--","Marker","*","MarkerSize",10)

quantiles = quantile(xtrace_errs(:,:,2),[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"LineStyle","-","Marker","x","MarkerSize",10)
quantiles = quantile(xtrace_iso_errs(:,:,2),[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),darkblue,"LineStyle","--","Marker","x","MarkerSize",10)

quantiles = quantile(xnystrace_errs(:,:,2),[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"LineStyle","-","Marker","o","MarkerSize",10)
quantiles = quantile(xnystrace_iso_errs(:,:,2),[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),darkorange,"LineStyle","--","Marker","o","MarkerSize",10)


set(gca,"YScale","log")

xlabel("Number of matvecs $s$")
ylabel("Relative error")
legend({"","Hutch++","","","","XTrace","","","","XNysTrace","",""},"Location","SouthWest")

exportgraphics(gcf,"../figs/fig_14_3_b.png")
saveas(gcf,"../figs/fig_14_3_b.fig")