clear;
clc;
close all;

a = dir(fullfile(pwd,'_*'));
for j = 1:length(a)
	foldername = a(j).name;
    disp(foldername);
    filename = ['./',foldername,'/',foldername,'_v.csv'];
    v = csvread(filename);
    filename = ['./',foldername,'/',foldername,'_v_x.csv'];
    vx = csvread(filename);
    filename = ['./',foldername,'/',foldername,'_v_y.csv'];
    vy = csvread(filename);
    filename = ['./',foldername,'/',foldername,'_v_z.csv'];
    vz = csvread(filename);
    dlmwrite(['./',foldername,'/',foldername,'_ACF.csv'], [], 'delimiter', ',');
    dlmwrite(['./',foldername,'/',foldername,'_dtheta.csv'], [], 'delimiter', ',');
    dlmwrite(['./',foldername,'/',foldername,'_dtheta_histogram.csv'], [], 'delimiter', ',');
    for i = 2:size(v,1)
        ACF =  sum(vx(i,:).*vx(1,:) + vy(i,:).*vy(1,:) + vz(i,:).*vz(1,:))/size(v,2);
        dtheta = 180*abs(acos((vx(i,:).*vx(i - 1,:) + vy(i,:).*vy(i - 1,:) + vz(i,:).*vz(i - 1,:))./(v(i,:) + eps)./(v(i - 1,:) + eps)*3600*3600))/pi;
    	dlmwrite(['./',foldername,'/',foldername,'_dtheta.csv'], dtheta, 'delimiter', ',','-append');
        dlmwrite(['./',foldername,'/',foldername,'_ACF.csv'], ACF', 'delimiter', ',','-append');
        dtheta_count = histcounts(dtheta,0:10:180);
        dlmwrite(['./',foldername,'/',foldername,'_dtheta_histogram.csv'], dtheta_count, 'delimiter', ',','-append');


    end
end