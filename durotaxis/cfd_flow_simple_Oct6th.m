
 clear;
 clc;
 close all;
mkdir('figure');
rmdir('figure','s');
mkdir('figure');


l = 1;
for kk = 1:5
   
 for k = 0:7
     for j = 1:3
         close all;
         filename = '';
         if kk == 1
             ECM_stiffness_gradient = 0;
             filename = [filename,'_dE00'];
         elseif kk == 2
             ECM_stiffness_gradient = 0.5;
             filename = [filename,'_dE05'];
         elseif kk == 3
             ECM_stiffness_gradient = 1;
             filename = [filename,'_dE10'];
         elseif kk == 4
             ECM_stiffness_gradient = 1.5;
             filename = [filename,'_dE15'];
         elseif kk == 5
             ECM_stiffness_gradient = 2;
             filename = [filename,'_dE20'];
         end
         if k == 0
            filename =[filename,'_00um_',num2str(j)];
         elseif k == 1
            filename = [filename,'_03um_',num2str(j)];
         elseif k == 2
             filename = [filename,'_05um_',num2str(j)];
         elseif k == 3
             filename = [filename,'_10um_',num2str(j)];
         elseif k == 4
             filename = [filename,'_15um_',num2str(j)];
         elseif k == 5
             filename = [filename,'_20um_',num2str(j)];
         elseif k == 6
             filename = [filename,'_25um_',num2str(j)];
         elseif k == 7
             filename = [filename,'_30um_',num2str(j)];
         end
         clearvars p x cross_pairs;
         p = cfd_initiation_flow_May17th_02mg(filename,k,kk);
         [p,x,cross_pairs,h] = cfd_FiberGenerator_flow_Sep26th(p,kk);
         if round(p.CollagenDensity*100)/100 ~= p.rho
%              fprintf('regenerate fibernetwork\n');
             p.n = floor(p.n/p.CollagenDensity*p.rho);
             [p,x,cross_pairs,h] = cfd_FiberGenerator_flow_Sep26th(p,kk);
         end
         % Here we design a gradient of E
        x.z_seg_mid = (x.z_seg_start + x.z_seg_end)/2;
        if ECM_stiffness_gradient == 0
            p.E = ones(size(x.z_seg_mid))*p.E;
        else
            p.E = (x.z_seg_mid - p.bottom_scale + p.domain_scale_z/2)/p.domain_scale_z*ECM_stiffness_gradient/2*p.E + p.E;
        end
    %     p.D = (p.D - 0.05*p.D) + (x.z_seg_mid - p.bottom_scale)*0.1*p.D/p.domain_scale_z;
        p.seg_diameter = p.seg_diameter;
        p.A = pi*p.seg_diameter.*p.seg_diameter/4;
        p.Irot = 1/64*pi*p.seg_diameter.*p.seg_diameter.*p.seg_diameter.*p.seg_diameter; % Second moment of inertia for circular cross section for all cartesian axes
        p.kbend = p.E.*p.Irot/p.fiber_length_mean;
         x_cell = cell_initiation_flow_May17th(p,kk);
         i = 1;
         h = 0;
        %  h = cfd_output_flow_May18th2017(p,x,x_cell,1,h);
        fprintf(['collagen density: ',num2str(p.CollagenDensity),'\n']);
         fprintf(['collagen fiber number: ',num2str(p.n),'\n']);
         fprintf(['ECM stiffness gradient percentage: ',num2str(ECM_stiffness_gradient*100),'\n']);
         fprintf(['flowspeed: ',num2str(k),'\n']);
         fprintf(['trial: ',num2str(j),'\n']);
        %  x.rotcenter = p.rotcenter;
         F = TractionForceCalculator(p,x,x_cell);
         F = ProtursionForceCalculator(F,p,x,x_cell);
         MSD = zeros();
         speed_magnitude = zeros();
         cell_theta = zeros();
         cell_phi = zeros();
         cell_positions_x = zeros(p.tnind + 1,x_cell.N);
         cell_positions_y = zeros(p.tnind + 1,x_cell.N);
         cell_positions_z = zeros(p.tnind + 1,x_cell.N);
         drag = 6*pi*repmat(x_cell.R,1,3)*p.viscosity.*(repmat(p.flow_v,x_cell.N,1)*1e-6 - x_cell.v);
        %  figure;
        % % axis hold;
        % plot3(x_cell.location_x,x_cell.location_y,x_cell.location_z,'r.','Markersize',12,'LineWidth',100);
        % axis([-20 20 -20 20 -20 20]);
        dlmwrite(['./',p.filename,'/',p.filename,'_x.csv'], x_cell.location_x', 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_y.csv'], x_cell.location_y', 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_z.csv'], x_cell.location_z', 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_v_x.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_v_y.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_v_z.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_MSD.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_v.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_theta.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_phi.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_v_theta.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_v_phi.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_T_theta.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_T_phi.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_T.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_T_fiberN.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_R_fiberN.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_Tx.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_Ty.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_Tz.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_Rz.csv'], [], 'delimiter', ',');
        dlmwrite(['./',p.filename,'/',p.filename,'_TRz.csv'], [], 'delimiter', ',');


        % hold on;
        for i = 1:p.tnind
                   %x_cell.v = zeros();
            for m = 1:100
                for n = 1:3
                    if abs(drag(m,n) + F.protrusion(m,n) + F.traction(m,n)) <= abs(F.frictionECM(m,n))
                       x_cell.v (m,n) = 0;
                    else
                    x_cell.v (m,n) = (drag (m,n) + F.protrusion (m,n) + F.traction (m,n) + F.frictionECM(m,n))./(p.viscosity_ECM*6*pi*7.5);
                   %x_cell.v (m,n) = (drag (m,n) + F.protrusion (m,n) + F.frictionECM(m,n))./(p.viscosity_ECM*6*pi*7.5);
                   %x_cell.v (m,n) = (drag (m,n) + F.protrusion (m,n) + F.traction (m,n))./(p.viscosity_ECM*6*pi*7.5);
                    end
            end
            end
            original_location_x = x_cell.location_x;
            original_location_y = x_cell.location_y;
            original_location_z = x_cell.location_z;

            cell_positions_x(i,:) = original_location_x';
            cell_positions_y(i,:) = original_location_y';
            cell_positions_z(i,:) = original_location_z';

            x_cell.location_x = x_cell.location_x + x_cell.v(:,1)*p.dt;
            x_cell.location_y = x_cell.location_y + x_cell.v(:,2)*p.dt;
            x_cell.location_z = x_cell.location_z + x_cell.v(:,3)*p.dt;

            MSD =  (sum((x_cell.location_x - x_cell.location_x_source).^2) + sum((x_cell.location_y - x_cell.location_y_source).^2) + sum((x_cell.location_z - x_cell.location_z_source).^2))/x_cell.N;

            speed_magnitude = sqrt((x_cell.v(:,1)).^2 + (x_cell.v(:,2)).^2 + (x_cell.v(:,3)).^2);

            F = TractionForceCalculator(p,x,x_cell);
            F = ProtursionForceCalculator(F,p,x,x_cell);

            v_theta = acos(x_cell.v(:,3)*p.dt./sqrt((x_cell.v(:,1)*p.dt).^2 + (x_cell.v(:,2)*p.dt).^2 + (x_cell.v(:,3)*p.dt).^2))*180/pi;
            v_phi = atan2(x_cell.v(:,2)*p.dt,x_cell.v(:,1)*p.dt)*180/pi;
            
            tract_theta = acos(F.traction(:,3)./sqrt((F.traction(:,1)).^2 + (F.traction(:,2)).^2 + (F.traction(:,3)).^2))*180/pi;
            tract_phi = atan2(F.traction(:,2),F.traction(:,1))*180/pi;
            
            cell_theta = acos(x_cell.location_z./sqrt((x_cell.location_x).^2 + (x_cell.location_y).^2 + (x_cell.location_z).^2))*180/pi;
            cell_phi = atan2(x_cell.location_y,x_cell.location_x)*180/pi;

        
            dlmwrite(['./',p.filename,'/',p.filename,'_x.csv'], x_cell.location_x','delimiter',',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_y.csv'], x_cell.location_y','delimiter',',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_z.csv'], x_cell.location_z','delimiter',',','-append');
            
            dlmwrite(['./',p.filename,'/',p.filename,'_T.csv'], sqrt(F.traction(:,1).^2+F.traction(:,2).^2+F.traction(:,3).^2)'*1e12,'delimiter',',','-append');

            dlmwrite(['./',p.filename,'/',p.filename,'_v_x.csv'], x_cell.v(:,1)','delimiter',',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_v_y.csv'], x_cell.v(:,2)','delimiter',',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_v_z.csv'], x_cell.v(:,3)','delimiter',',','-append');

            dlmwrite(['./',p.filename,'/',p.filename,'_MSD.csv'], MSD','delimiter',',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_v.csv'], 3600*speed_magnitude','delimiter',',','-append');

            dlmwrite(['./',p.filename,'/',p.filename,'_theta.csv'],cell_theta','delimiter',',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_phi.csv'],cell_phi','delimiter',',','-append');
            
            dlmwrite(['./',p.filename,'/',p.filename,'_T_theta.csv'],tract_theta','delimiter',',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_T_phi.csv'],tract_phi','delimiter',',','-append');
            
         
            
            dlmwrite(['./',p.filename,'/',p.filename,'_v_theta.csv'], v_theta','delimiter',',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_v_phi.csv'],v_phi','delimiter',',','-append');
            
            dlmwrite(['./',p.filename,'/',p.filename,'_T_fiberN.csv'], F.T_fiberN,'delimiter',',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_R_fiberN.csv'], F.R_fiberN,'delimiter',',','-append');

            
            dlmwrite(['./',p.filename,'/',p.filename,'_Tx.csv'], F.traction(:,1)', 'delimiter', ',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_Ty.csv'], F.traction(:,2)', 'delimiter', ',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_Tz.csv'], F.traction(:,3)', 'delimiter', ',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_Rz.csv'], F.frictionECM(:,3)', 'delimiter', ',','-append');
            dlmwrite(['./',p.filename,'/',p.filename,'_TRz.csv'], F.traction(:,3)' + F.frictionECM(:,3)', 'delimiter', ',','-append');


        end
        
        
        
        
        
        
        
        l = l + 1;
        
     end
     
 end
 
 end

