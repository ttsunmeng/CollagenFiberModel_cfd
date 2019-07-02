% Having problem calculating the strains.
% Please add the reorientation


function [x,h] = cfd_output_Aug23th(p,x,F,findex,cross_pairs,i,h)


if i == 1
    dlmwrite(['./',p.filename,'/',p.filename,'move_border.csv'], [], 'delimiter', ',');
    dlmwrite(['./',p.filename,'/',p.filename,'z_Stress.csv'], [], 'delimiter', ',');
    dlmwrite(['./',p.filename,'/',p.filename,'crosslink.csv'], [i,length(cross_pairs)], 'delimiter', ',');
    dlmwrite(['./',p.filename,'/',p.filename,'_energy.csv'], [], 'delimiter', ',');
    dlmwrite(['./',p.filename,'/',p.filename,'_num.csv'], [], 'delimiter', ',');

    
    
    if p.vis_in_VMD == 1
        vmd_psf_generator(findex.fiber_nodes_num,p.filename);
    	vmd_pdb_generator(1,findex.fiber_nodes_num,x.x_seg_cor_glo,x.y_seg_cor_glo,x.z_seg_cor_glo,p.filename);
    end
    
    x.findex_x_org = x.x_seg_cor_glo(findex.findex_right_nind) - x.x_seg_cor_glo(findex.findex_left_nind);
    x.findex_y_org = x.y_seg_cor_glo(findex.findex_right_nind) - x.y_seg_cor_glo(findex.findex_left_nind);
    x.findex_z_org = x.z_seg_cor_glo(findex.findex_right_nind) - x.z_seg_cor_glo(findex.findex_left_nind);

    findex_length_tmp = sqrt((x.findex_x_org).^2 + (x.findex_y_org).^2 + (x.findex_z_org).^2);
    theta_tmp = acos(x.findex_z_org./(findex_length_tmp + eps))*180/pi; % 0 to 180
    phi_tmp = atan(x.findex_y_org./(x.findex_x_org + eps))*180/pi; % 0 to 180;
    z_middle = (x.z_seg_cor_glo(findex.findex_right_nind) + x.z_seg_cor_glo(findex.findex_left_nind))/2;
    thetaCount_all = zeros(19*20,1);
    phiCount_all = zeros(19*20,1);
    for i = 1:20
        [thetaCount,~] = hist(theta_tmp(z_middle > 300 - i*5 & z_middle <= 300 - (i - 1)*5),0:10:180);
        thetaCount_all((i-1)*19 + 1:i*19) = thetaCount';
        [phiCount,~] = hist(phi_tmp(z_middle > 300 - i*5 & z_middle <= 300 - (i - 1)*5),0:10:180);
        phiCount_all((i-1)*19 + 1:i*19) = phiCount';
        
    end
        [lengthCount,~] = hist(findex_length_tmp,0:0.05:5);
        [DensityCount,~] = hist(z_middle,200:0.5:300);

    dlmwrite(['./',p.filename,'/',p.filename,'orient_theta.csv'], [i thetaCount_all'], 'delimiter', ',');
    dlmwrite(['./',p.filename,'/',p.filename,'orient_phi.csv'], [i phiCount_all'], 'delimiter', ',');
    dlmwrite(['./',p.filename,'/',p.filename,'length.csv'], [i lengthCount], 'delimiter', ',');
    dlmwrite(['./',p.filename,'/',p.filename,'density.csv'], [i DensityCount], 'delimiter', ',');


    if p.vis_in_matlab == 1 || p.safe_vis == 1  
        x_vis_start = x.x_seg_cor_glo;
        x_vis_end = x.x_seg_cor_glo;
        x_vis_start(findex.glo_right_ind) = NaN;
        x_vis_end(findex.glo_left_ind) = NaN;
        x_vis_start(isnan(x_vis_start)) = [];
        x_vis_end(isnan(x_vis_end)) = [];

        y_vis_start = x.y_seg_cor_glo;
        y_vis_end = x.y_seg_cor_glo;
        y_vis_start(findex.glo_right_ind) = NaN;
        y_vis_end(findex.glo_left_ind) = NaN;
        y_vis_start(isnan(y_vis_start)) = [];
        y_vis_end(isnan(y_vis_end)) = [];

        z_vis_start = x.z_seg_cor_glo;
        z_vis_end = x.z_seg_cor_glo;
        z_vis_start(findex.glo_right_ind) = NaN;
        z_vis_end(findex.glo_left_ind) = NaN;
        z_vis_start(isnan(z_vis_start)) = [];
        z_vis_end(isnan(z_vis_end)) = [];
        %% Make it adjustable with the disp!!
        if p.vis_in_matlab == 1
            ff = figure('Position',[0 0 750 800]);
        else
            ff = figure('Visible','off','Position',[0 0 750 800]);
        end
        h = plot3([transpose(x_vis_start);transpose(x_vis_end)],[transpose(y_vis_start);transpose(y_vis_end)],[transpose(z_vis_start);transpose(z_vis_end)],'b-');
        xlim([p.bottom_scale - 5,p.bottom_scale + p.cube_scale + 5]);
        ylim([p.bottom_scale - 5,p.bottom_scale + p.cube_scale + 5]);
        zlim([p.bottom_scale - 5,p.bottom_scale + p.cube_scale + 20]);
        
        
        if p.safe_vis == 1
            saveas(ff,['./',p.filename,'/final_fig_',num2str(i)],'png');
        end
    else
        h = 0;
    end
    
elseif i >= 2
    dlmwrite(['./',p.filename,'/',p.filename,'move_border.csv'], [i,p.move_border(i-1)], 'delimiter', ',','-append');
    dlmwrite(['./',p.filename,'/',p.filename,'z_Stress.csv'], [i,F.stress_z], 'delimiter', ',','-append');
    dlmwrite(['./',p.filename,'/',p.filename,'crosslink.csv'], [i,size(cross_pairs,1)], 'delimiter', ',','-append');
    dlmwrite(['./',p.filename,'/',p.filename,'_energy.csv'], [i,F.Estretch(i),F.Ebend(i),F.Ecrosslink(i),F.Estretch(i)+F.Ebend(i)+F.Ecrosslink(i)], 'delimiter', ',','-append');
    dlmwrite(['./',p.filename,'/',p.filename,'_num.csv'], [i,size(findex.interface,1),sum(F.F_stretch_z>0&x.z_seg_cor_glo >= p.move_border(i) - p.interface_depth),sum(F.F_stretch_z<0&x.z_seg_cor_glo >= p.move_border(i) - p.interface_depth),sum(F.F_bend_z>0&x.z_seg_cor_glo >= p.move_border(i) - p.interface_depth),sum(F.F_bend_z<0&x.z_seg_cor_glo >= p.move_border(i) - p.interface_depth),F.stress_z/(p.cube_scale^2)*1e9], 'delimiter', ',','-append');
    x.findex_x_org = x.x_seg_cor_glo(findex.findex_right_nind) - x.x_seg_cor_glo(findex.findex_left_nind);
    x.findex_y_org = x.y_seg_cor_glo(findex.findex_right_nind) - x.y_seg_cor_glo(findex.findex_left_nind);
    x.findex_z_org = x.z_seg_cor_glo(findex.findex_right_nind) - x.z_seg_cor_glo(findex.findex_left_nind);

    findex_length_tmp = sqrt((x.findex_x_org).^2 + (x.findex_y_org).^2 + (x.findex_z_org).^2);
    theta_tmp = acos(x.findex_z_org./(findex_length_tmp + eps))*180/pi;
    phi_tmp = atan(x.findex_y_org./(x.findex_x_org + eps))*180/pi; % 0 to 180;
    z_middle = (x.z_seg_cor_glo(findex.findex_right_nind) + x.z_seg_cor_glo(findex.findex_left_nind))/2;
    thetaCount_all = zeros(19*20,1);
    phiCount_all = zeros(19*20,1);
    for i = 1:20
        [thetaCount,~] = hist(theta_tmp(z_middle > 300 - i*5 & z_middle <= 300 - (i - 1)*5),0:10:180);
        thetaCount_all((i-1)*19 + 1:i*19) = thetaCount';
        [phiCount,~] = hist(phi_tmp(z_middle > 300 - i*5 & z_middle <= 300 - (i - 1)*5),0:10:180);
        phiCount_all((i-1)*19 + 1:i*19) = phiCount';
    end
        [lengthCount,~] = hist(findex_length_tmp,0:0.05:5);
        [DensityCount,~] = hist(z_middle,200:0.5:300);

    dlmwrite(['./',p.filename,'/',p.filename,'orient_theta.csv'], [i thetaCount_all'], 'delimiter', ',','-append');
    dlmwrite(['./',p.filename,'/',p.filename,'orient_phi.csv'], [i phiCount_all'], 'delimiter', ',','-append');
    dlmwrite(['./',p.filename,'/',p.filename,'length.csv'], [i lengthCount], 'delimiter', ',','-append');
    dlmwrite(['./',p.filename,'/',p.filename,'density.csv'], [i DensityCount], 'delimiter', ',','-append');
     
    if p.vis_in_matlab == 1 || p.safe_vis == 1
        x_vis_start = x.x_seg_cor_glo;
        x_vis_end = x.x_seg_cor_glo;
        x_vis_start(findex.glo_right_ind) = NaN;
        x_vis_end(findex.glo_left_ind) = NaN;
        x_vis_start(isnan(x_vis_start)) = [];
        x_vis_end(isnan(x_vis_end)) = [];

        y_vis_start = x.y_seg_cor_glo;
        y_vis_end = x.y_seg_cor_glo;
        y_vis_start(findex.glo_right_ind) = NaN;
        y_vis_end(findex.glo_left_ind) = NaN;
        y_vis_start(isnan(y_vis_start)) = [];
        y_vis_end(isnan(y_vis_end)) = [];

        z_vis_start = x.z_seg_cor_glo;
        z_vis_end = x.z_seg_cor_glo;
        z_vis_start(findex.glo_right_ind) = NaN;
        z_vis_end(findex.glo_left_ind) = NaN;
        z_vis_start(isnan(z_vis_start)) = [];
        z_vis_end(isnan(z_vis_end)) = [];
        if mod(i,200) == 0
            close all;
            if p.vis_in_matlab == 1
                ff = figure('Position',[0 0 750 800]);
            else
                
                ff = figure('Visible','off','Position',[0 0 750 800]);
            end
            h = plot3([transpose(x_vis_start);transpose(x_vis_end)],[transpose(y_vis_start);transpose(y_vis_end)],[transpose(z_vis_start);transpose(z_vis_end)],'b-');
            xlim([p.bottom_scale - 5,p.bottom_scale + p.cube_scale + 5]);
            ylim([p.bottom_scale - 5,p.bottom_scale + p.cube_scale + 5]);
            zlim([p.bottom_scale - 5,p.bottom_scale + p.cube_scale + 20]);
            
            if p.safe_vis == 1
                saveas(ff,['./',p.filename,'/final_fig_',num2str(i)],'png');
            end
        end
    end
    
    if p.vis_in_VMD == 1
%           vmd_pdb_generator(i,findex.fiber_nodes_num,x.x_seg_cor_glo,x.y_seg_cor_glo,x.z_seg_cor_glo,p.filename);

        if mod(i,100) == 0
            vmd_pdb_generator(i/100+1,findex.fiber_nodes_num,x.x_seg_cor_glo,x.y_seg_cor_glo,x.z_seg_cor_glo,p.filename);
        end
    else
        h = 0;
    end

    
end
    
