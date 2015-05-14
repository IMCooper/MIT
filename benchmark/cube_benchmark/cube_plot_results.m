%bench2:
% filenames = {'~/Desktop/temp_mnt/code_library/benchmark/cube_benchmark/results/bench2_diagaxis_ptfield_p0.out' ...
%     '~/Desktop/temp_mnt/code_library/benchmark/cube_benchmark/results/bench2_diagaxis_ptfield_p1.out' ...
%     '~/Desktop/temp_mnt/code_library/benchmark/cube_benchmark/results/bench2_diagaxis_ptfield_p2.out' ...
%     '~/Desktop/temp_mnt/code_library/benchmark/cube_benchmark/results/bench2_diagaxis_ptfield_p3.out' };
%bench3:
%  filenames = {'~/Desktop/temp_mnt/code_library/benchmark/cube_benchmark/results/bench3_diagaxis_ptfield_p0.out' ...
%      '~/Desktop/temp_mnt/code_library/benchmark/cube_benchmark/results/bench3_diagaxis_ptfield_p1.out' ...
%      '~/Desktop/temp_mnt/code_library/benchmark/cube_benchmark/results/bench3_diagaxis_ptfield_p2.out' ...
%      '~/Desktop/temp_mnt/code_library/benchmark/cube_benchmark/results/bench3_diagaxis_ptfield_p3.out'};

% temp:
filenames = {'~/Desktop/temp_mnt/new_code/benchmark/cube_benchmark/build/delete3_diagaxis_ptfield_p0.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/cube_benchmark/build/delete3_diagaxis_ptfield_p1.out' ...
    '~/Desktop/temp_mnt/new_code/benchmark/cube_benchmark/build/delete3_diagaxis_ptfield_p2.out'};% ...
    %'~/Desktop/temp_mnt/new_code/benchmark/cube_benchmark/build/delete2_diagaxis_ptfield_p2.out'};

line_colour = {'.-r' '.-k' '.-g' '.-m'};

% load the first:
data = load(filenames{1});
r = data(:,1);
Ee = data(:,8:10)' + 1i*data(:,11:13)';
He = data(:,20:22)' + 1i*data(:,23:25)';
He_pt = data(:,32:34)' + 1i*data(:,35:37)';

% REMOVED:
% Adjust for perturbed:
% H0=[0;0;1];
% for i=1:length(r)
%     He(:,i) = He(:,i) - H0;
% end

% Magnitudes:
magEe=zeros(length(r),1);
magHe=zeros(length(r),1);
magHe_pt=zeros(length(r),1);
for i=1:length(r)
    magEe(i) = norm(Ee(:,i));
    magHe(i) = norm(He(:,i));
    magHe_pt(i) = norm(He_pt(:,i));
end
h1=figure;
semilogy(r,magEe);
h2=figure;
semilogy(r,magHe);
h3=figure;
semilogy(r,magHe_pt);

magEa=zeros(length(r),length(filenames));
magHa=zeros(length(r),length(filenames));
magHa_pt=zeros(length(r),length(filenames));
for f=1:length(filenames)
    data=load(filenames{f});
    % approx H:
    Ea = data(:,2:4)' + 1i*data(:,5:7)';
    Ha = data(:,14:16)' + 1i*data(:,17:19)';
    Ha_pt = data(:,26:28)' + 1i*data(:,29:31)';
    % REMOVED:
    % Adjust for perturbed:
%     for i=1:length(r)
%         Ha(:,i) = Ha(:,i)-H0;
%     end
    
    for i=1:length(r)
        magEa(i,f) = norm(Ea(:,i));
        magHa(i,f) = norm(Ha(:,i));
        magHa_pt(i,f) = norm(Ha_pt(:,i));
    end
    % Add to plots
    figure(h1);
    hold on;
    semilogy(r,magEa(:,f),line_colour{f});
    hold off;
    figure(h2);
    hold on;
    semilogy(r,magHa(:,f),line_colour{f});
    hold off;
    figure(h3);
    hold on;
    semilogy(r,magHa_pt(:,f),line_colour{f});
    hold off;
end
figure(h1);
legend('PT', 'p=0', 'p=1', 'p=2', 'p=3', 'Location', 'Best');
figure(h2);
legend('PT', 'p=0', 'p=1', 'p=2', 'p=3', 'Location', 'Best');
figure(h3);
legend('PT', 'p=0', 'p=1', 'p=2', 'p=3', 'Location', 'Best');

% saveas(h1,'cube_direct_magE','fig');
% saveas(h1,'cube_direct_magE','epsc2');
% saveas(h2,'cube_direct_magH','fig');
% saveas(h2,'cube_direct_magH','epsc2');
% saveas(h3,'cube_direct_magHp','fig');
% saveas(h3,'cube_direct_magHp','epsc2');
