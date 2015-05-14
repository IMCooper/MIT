% Script to plot the older outputs
% The newer output files contain more info (both the full h field and the
% perturbed part only), where as we need to calculate the perturbed one
% manually with these older files.

filenames = {'results/sphere_bench134_old/bench1_diag_ptfield_p0.out' ...
    'results/sphere_bench134_old/bench3_diag_ptfield_p0.out' ...
    'results/sphere_bench134_old/bench4_diag_ptfield_p0.out'};

line_colour = {'r' 'k' 'g' 'm'};

% load the first:
data = load(filenames{1});
r = data(:,1);
He = data(:,20:22)' + 1i*data(:,23:25)';

% Adjust for perturbed:
H0=[0;0;1];
for i=1:length(r)
    He(:,i) = He(:,i)-H0;
end

% Magnitudes:
magHe=zeros(length(r),1);
magHe_pt=zeros(length(r),1);
for i=1:length(r)
    magHe(i) = norm(He(:,i));
end
figure;
semilogy(r,magHe);
hold on;

magHa=zeros(length(r),length(filenames));
magHa_pt=zeros(length(r),length(filenames));
for f=1:length(filenames)
    data=load(filenames{f});
    % approx H:
    Ha = data(:,14:16)' + 1i*data(:,17:19)';
    % Adjust for perturbed:
    for i=1:length(r)
        Ha(:,i) = Ha(:,i)-H0;
    end
    
    for i=1:length(r)
        magHa(i,f) = norm(Ha(:,i));
    end
    % Add to plot
    semilogy(r,magHa(:,f),line_colour{f});
end

legend('PT', 'M1', 'M2', 'M3', 'M4');