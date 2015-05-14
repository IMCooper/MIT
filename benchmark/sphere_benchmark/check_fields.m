function data = check_fields(input_file,savename)

plotcomps=0;

data=load(input_file);
r = data(:,1);
% approx sol:
Ea = data(:,2:4)' + 1i*data(:,5:7)';
Ha = data(:,14:16)' + 1i*data(:,17:19)';
Ha_pt = data(:,26:28)' + 1i*data(:,29:31)';

% exact sol:
Ee = data(:,8:10)' + 1i*data(:,11:13)';
He = data(:,20:22)' + 1i*data(:,23:25)';
He_pt = data(:,32:34)' + 1i*data(:,35:37)';

% Magnitudes:
magEe=zeros(length(r),1);
magHe=magEe;
magHe_pt=magEe;
magEa=magEe;
magHa=magEe;
magHa_pt=magEe;
for i=1:length(r)
    magEe(i) = norm(Ee(:,i));
    magEa(i) = norm(Ea(:,i));    
    magHe(i) = norm(He(:,i));
    magHa(i) = norm(Ha(:,i));
    magHe_pt(i) = norm(He_pt(:,i));
    magHa_pt(i) = norm(Ha_pt(:,i));
    
end


mu0=4*pi*1e-7;
omega=2*pi*50;clc

factor=1;

component = {'x' 'y' 'z'};

% close all;
if plotcomps
for dim_E=1:3
    figure;
    plot(r,real(Ee(dim_E,:))); hold on; plot(r,real(Ea(dim_E,:)),'r')
    
    if exist('savename','var')
        saveas(gcf,strcat(savename,'_reA',component{dim_E}),'fig');
        saveas(gcf,strcat(savename,'_reA',component{dim_E}),'epsc2');
    end
    
    figure;
    plot(r,imag(Ee(dim_E,:))); hold on; plot(r,imag(Ea(dim_E,:)),'r')
    if exist('savename','var')
        saveas(gcf,strcat(savename,'_imA',component{dim_E}),'fig');
        saveas(gcf,strcat(savename,'_imA',component{dim_E}),'epsc2');
    end
end
for dim_H=1:3
    figure;
    plot(r,real(He(dim_H,:))); hold on; plot(r,real(Ha(dim_H,:)),'r')
    if exist('savename','var')
        saveas(gcf,strcat(savename,'_recurlA',component{dim_H}),'fig');
        saveas(gcf,strcat(savename,'_recurlA',component{dim_H}),'epsc2');
    end
    figure;
    plot(r,imag(He(dim_H,:))); hold on; plot(r,imag(Ha(dim_H,:)),'r')
    if exist('savename','var')
        saveas(gcf,strcat(savename,'_imcurlA',component{dim_H}),'fig');
        saveas(gcf,strcat(savename,'_imcurlA',component{dim_H}),'epsc2');
    end
end
end

% Plot magnitudes:
figure;
semilogy(r, magHe, r, magHa);
figure;
semilogy(r,magHe_pt, r, magHa_pt);


end