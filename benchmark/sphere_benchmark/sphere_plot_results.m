function sphere_plot_results(filenames)

% filenames = {'results/sphere_bench1345_p0/bench1_diag_xaxis_ptfield_p0.out' ...
%     'results/sphere_bench1345_p0/bench3_diag_xaxis_ptfield_p0.out' ...
%     'results/sphere_bench1345_p0/bench4_diag_xaxis_ptfield_p0.out' ...
%     'results/sphere_bench1345_p0/bench5_diag_xaxis_ptfield_p0.out' };

% filenames = {'results/sphere_meshrefinement_run1/bench1_diag_ptfield_p0.out' ...
%     'results/sphere_meshrefinement_run1/bench3_diag_ptfield_p0.out' ...
%     'results/sphere_meshrefinement_run1/bench4_diag_ptfield_p0.out' ...
%     'results/sphere_meshrefinement_run1/bench5_diag_ptfield_p0.out' };

% filenames = {'results/bench1_xaxis_ptfield_p0.out' ...
%     'results/bench3_xaxis_ptfield_p0.out' ...
%     'results/bench4_xaxis_ptfield_p0.out'};


line_colour = {'.-r' '+-k' '*-g' '^-m' 'd-c' 'x-y'};

% load the first:
data = load(filenames{1});
r = data(:,1);
He = data(:,20:22)' + 1i*data(:,23:25)';
He_pt = data(:,32:34)' + 1i*data(:,35:37)';

% Adjust for perturbed:
% H0=[0;0;1];
% for i=1:length(r)
%     He(:,i) = He(:,i);
% end

% Magnitudes:
magHe=zeros(length(r),1);
magHe_pt=zeros(length(r),1);
for i=1:length(r)
    magHe(i) = norm(He(:,i));
    magHe_pt(i) = norm(He_pt(:,i));
end
h1=figure;
semilogy(r,magHe,'LineWidth',2);
h2=figure;
semilogy(r,magHe_pt,'LineWidth',2);


magHa=zeros(length(r),length(filenames));
magHa_pt=zeros(length(r),length(filenames));
for f=1:length(filenames)
    data=load(filenames{f});
    % approx H:
    Ha = data(:,14:16)' + 1i*data(:,17:19)';
    Ha_pt = data(:,26:28)' + 1i*data(:,29:31)';
    % Adjust for perturbed:
    for i=1:length(r)
        Ha(:,i) = Ha(:,i);
    end
    
    for i=1:length(r)
        magHa(i,f) = norm(Ha(:,i));
        magHa_pt(i,f) = norm(Ha_pt(:,i));
    end
    % Add to plot
    figure(h1);
    hold on;
    semilogy(r,magHa(:,f),line_colour{f});
    hold off;
    figure(h2);
    hold on;
    semilogy(r,magHa_pt(:,f),line_colour{f});
    hold off;
    
end


for i=[h1 h2]
    figure(i);
    
% mesh ref:
% legend('PT', 'M1', 'M2', 'M3', 'M4');

% p-ref:
% legend('PT', 'p=0', 'p=1', 'p=2', 'p=3');

legend('PT', 'p=1, g=1', 'p=2, g=1', 'p=1, g=2', 'p=2, g=2');
end


end