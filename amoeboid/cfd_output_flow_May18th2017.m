% Having problem calculating the strains.
% Please add the reorientation


function [x,h] = cfd_output_flow_May18th2017(p,x,x_cell,i,h)

    close all;
    if p.vis_in_matlab == 1
        ff = figure('Position',[0 0 750 750*8/4]);
    else
        ff = figure('Visible','off','Position',[0 0 750 750*6/4]);
    end
    axis([p.bottom_scale - 5,p.bottom_scale + p.domain_scale_y + 5,p.bottom_scale - 5,p.bottom_scale + p.domain_scale_y + 5,p.bottom_scale - 5,p.bottom_scale + p.domain_scale_z + 20]);
    [xsph,ysph,zsph]=sphere(35);
    crad = 7.5;
    h = plot3([transpose(x.x_seg_start);transpose(x.x_seg_end)],[transpose(x.y_seg_start);transpose(x.y_seg_end)],[transpose(x.z_seg_start);transpose(x.z_seg_end)],'b-','linewidth',2);
    h.Color(4)=0.1;
    
%     camtarget([0 0 0]);
%     zoom(1);
%     light('Position',[-10,-10,0],'Style','infinite');
    hold on;
    surf(x_cell.location_x+crad*xsph,x_cell.location_y+crad*ysph,x_cell.location_z+crad*zsph,'facecolor',[1 0 0],'edgecolor','none');        
    if p.safe_vis == 1
        saveas(ff,['./',p.filename,'/final_fig_',num2str(i)],'png');
    end
 

    
    

    
