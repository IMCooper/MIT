close all;
saveflag = 0;
highermapping1=2;
highermapping2=3;
for i=1:2
%     files{i}=strcat('newbuild/newelement_direct/cube_diagaxis_ptfield_p',num2str(i),'.out');
    files{i}=strcat('newbuild/mapping',num2str(highermapping1),'/cube_diagaxis_ptfield_p',num2str(i),'.out');
    files{i+2}=strcat('newbuild/mapping',num2str(highermapping2),'/cube_diagaxis_ptfield_p',num2str(i),'.out');
end
sphere_plot_results(files)
if saveflag
    figure(2);
    saveas(gcf,'sphere_low_v_quad_diag','fig');
    saveas(gcf,'sphere_low_v_quad_diag','epsc2');
    close all;
end


for i=1:2
%     files{i}=strcat('newbuild/newelement_direct/cube_xaxis_ptfield_p',num2str(i),'.out');
    files{i}=strcat('newbuild/mapping',num2str(highermapping1),'/cube_xaxis_ptfield_p',num2str(i),'.out');
    files{i+2}=strcat('newbuild/mapping',num2str(highermapping2),'/cube_xaxis_ptfield_p',num2str(i),'.out');
end
sphere_plot_results(files)
if saveflag
    figure(2);
    saveas(gcf,'sphere_low_v_quad_xaxis','fig');
    saveas(gcf,'sphere_low_v_quad_xaxis','epsc2');
    close all;
end

for i=1:2
%     files{i}=strcat('newbuild/newelement_direct/cube_yaxis_ptfield_p',num2str(i),'.out');
    files{i}=strcat('newbuild/mapping',num2str(highermapping1),'/cube_yaxis_ptfield_p',num2str(i),'.out');
    files{i+2}=strcat('newbuild/mapping',num2str(highermapping2),'/cube_yaxis_ptfield_p',num2str(i),'.out');
end
sphere_plot_results(files)
if saveflag
    figure(2);
    saveas(gcf,'sphere_low_v_quad_yaxis','fig');
    saveas(gcf,'sphere_low_v_quad_yaxis','epsc2');
    close all;
end

for i=1:2
%     files{i}=strcat('newbuild/newelement_direct/cube_zaxis_ptfield_p',num2str(i),'.out');
    files{i}=strcat('newbuild/mapping',num2str(highermapping1),'/cube_zaxis_ptfield_p',num2str(i),'.out');
    files{i+2}=strcat('newbuild/mapping',num2str(highermapping2),'/cube_zaxis_ptfield_p',num2str(i),'.out');
end
sphere_plot_results(files)
if saveflag
    figure(2);
    saveas(gcf,'sphere_low_v_quad_zaxis','fig');
    saveas(gcf,'sphere_low_v_quad_zaxis','epsc2');
    close all;
end

clear files;
for i=1:3
    files{i}=strcat('newbuild/mapping',num2str(i),'/cube_zaxis_ptfield_p',num2str(1),'.out');
end
sphere_plot_results(files)
if saveflag
    figure(2);
    saveas(gcf,'sphere_low_v_quad_zaxis','fig');
    saveas(gcf,'sphere_low_v_quad_zaxis','epsc2');
    close all;
end


clear files;
for i=1:2
    files{i}=strcat('newbuild/mapping1/cube_yaxis_ptfield_p',num2str(i),'.out');
    files{i+2}=strcat('newbuild/mappingstaticq1/cube_yaxis_ptfield_p',num2str(i),'.out');
end
sphere_plot_results(files)
if saveflag
    figure(2);
    saveas(gcf,'sphere_low_v_quad_zaxis','fig');
    saveas(gcf,'sphere_low_v_quad_zaxis','epsc2');
    close all;
end

close all;
clear files;
axis='diag';
for i=1:2
    files{i}=strcat('newbuild/cube_',axis,'axis_ptfield_p',num2str(i),'.out');
    files{i+2}=strcat('newbuild/mapping2/cube_',axis,'axis_ptfield_p',num2str(i),'.out');
end
sphere_plot_results(files)

close all;
clear files;
axis='diag';
p=4;
for g=1:6
    files{g}=strcat('~/Desktop/temp_mnt/newbuild/mapping',num2str(g),'/sphere_',axis,'axis_ptfield_p',num2str(p),'.out');
end
sphere_plot_results(files)
close(1);
legend('PT', 'g=1', 'g=2', 'g=3', 'g=4', 'g=5', 'g=6');

close all;
clear files;
axis='diag';
g=2;
for p=0:4
    files{p+1}=strcat('~/Desktop/temp_mnt/newbuild/mapping',num2str(g),'/sphere_',axis,'axis_ptfield_p',num2str(p),'.out');
end
sphere_plot_results(files)
close(1);
legend('PT', 'g=1', 'g=2', 'g=3', 'g=4', 'g=5', 'g=6');