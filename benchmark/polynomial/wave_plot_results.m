% temp:
filenames = {'~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark/results/sphere_check/sphere1_fields_0.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark/results/sphere_check/sphere2_fields_0.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark/results/sphere_check/sphere3_fields_0.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark/results/sphere_check/sphere4_fields_0.out' };
% Temp 2:
filenames = {'~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/cube_benchmark2/cube_fields_p0.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/cube_benchmark2/cube_fields_p1.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/cube_benchmark2/cube_fields_p2.out'};% ...
    %'~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/cube_benchmark2/cube1_fields_p3.out' ...
    %'~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/cube_benchmark2/cube1_fields_p4.out' };

% Temp 3:
% filenames = {'~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/deal_cube/cube64/cube1_fields_p0.out' ...
% '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/deal_cube/cube64/cube1_fields_p1.out' ...
% '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/deal_cube/cube64/cube1_fields_p2.out' ...
% '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/deal_cube/cube64/cube1_fields_p3.out'};

% Temp 2:
filenames = {'~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/deal_cube/cube_same_as_cubit/cube_fields_p0.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/deal_cube/cube_same_as_cubit/cube_fields_p1.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/deal_cube/cube_same_as_cubit/cube_fields_p2.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/waveprop_benchmark2/results/deal_cube/cube_same_as_cubit/cube_fields_p3.out'};



line_colour = {'.-r' '.-k' '.-g' '.-m' '.-y'};

% load the first:
data = load(filenames{1});
r = data(:,1);
Ee = data(:,8:10)' + 1i*data(:,11:13)';
He = data(:,20:22)' + 1i*data(:,23:25)';

% Magnitudes:
magHe=zeros(length(r),1);
magEe=zeros(length(r),1);
for i=1:length(r)
    magEe(i) = norm(Ee(:,i));
    magHe(i) = norm(He(:,i));    
end
h1=figure;
plot(r,magEe);
h2=figure;
plot(r,magHe);


magEe=zeros(length(r),length(filenames));
magHa=zeros(length(r),length(filenames));
for f=1:length(filenames)
    data=load(filenames{f});
    % approx H:
    Ea = data(:,2:4)' + 1i*data(:,5:7)';
    Ha = data(:,14:16)' + 1i*data(:,17:19)';
    
    for i=1:length(r)
        magEa(i,f) = norm(Ea(:,i));
        magHa(i,f) = norm(Ha(:,i));
    end
    % Add to plots
    figure(h1);
    hold on;
    plot(r,magEa(:,f),line_colour{f});
    hold off;
    figure(h2);
    hold on;
    plot(r,magHa(:,f),line_colour{f});
    hold off;
end
% figure(h1);
% legend('exact', 'M1', 'M2', 'M3', 'M4', 'Location', 'Best');
% figure(h2);
% legend('exact', 'M1', 'M2', 'M3', 'M4', 'Location', 'Best');

figure(h1);
legend('exact', 'p=0', 'p=1', 'p=2', 'p=3', 'Location', 'Best');
figure(h2);
legend('exact', 'p=0', 'p=1', 'p=2', 'p=3', 'Location', 'Best');
% 
saveas(h1,'dealcube_benchmark2_magE','fig');
saveas(h1,'dealcube_benchmark2_magE','epsc2');
saveas(h2,'dealcube_benchmark2_magH','fig');
saveas(h2,'dealcube_benchmark2_magH','epsc2');