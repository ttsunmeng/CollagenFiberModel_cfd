
 clear;
 clc;
 close all;
mkdir('figure');
rmdir('figure','s');
mkdir('figure');

kk = 2;

l = 1;

 for k = 1:1 %[0,1,7]
     
         
	
    for j = 1:1


        for CCL21_block = 0:0
            filename = 'cfd_flow_experiment_';
            if j == 1
             if CCL21_block == 0
                if k == 0
                    dlmwrite('cfd_00um_streamline_control.csv', [],'delimiter',',');
                    dlmwrite('cfd_00um_directional_control.csv', [],'delimiter',',');
                elseif k == 1
                    dlmwrite('cfd_03um_streamline_control.csv', [],'delimiter',',');
                    dlmwrite('cfd_03um_directional_control.csv', [],'delimiter',',');
                elseif k == 2
                    dlmwrite('cfd_3um_streamline_control.csv', [],'delimiter',',');
                    dlmwrite('cfd_3um_directional_control.csv', [],'delimiter',',');
                end
            elseif CCL21_block == 1
                if k == 0
                    dlmwrite('cfd_00um_streamline_CCL21block.csv', [],'delimiter',',');
                    dlmwrite('cfd_00um_directional_CCL21block.csv', [],'delimiter',',');
                elseif k == 1
                    dlmwrite('cfd_03um_streamline_CCL21block.csv', [],'delimiter',',');
                    dlmwrite('cfd_03um_directional_CCL21block.csv', [],'delimiter',',');
                elseif k == 2
                    dlmwrite('cfd_3um_streamline_CCL21block.csv', [],'delimiter',',');
                    dlmwrite('cfd_3um_directional_CCL21block.csv', [],'delimiter',',');
                end                   

             end
            end
             close all;
             
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
             if CCL21_block == 0
                filename =[filename,'_control'];
             elseif CCL21_block == 1
                filename = [filename,'_CCL21block'];
             end
             clearvars p x cross_pairs;
             p = cfd_initiation_flow_May17th_02mg(filename,k,kk);
             p.CCL21_block = CCL21_block;
             if CCL21_block == 1
             	p.temp_pool = (0:180)';
             end
             fprintf(['flowspeed: ',num2str(p.flow_v_z),'\n']);
             fprintf(['trial: ',num2str(j),'\n']);
             fprintf(['have CCL21 black? ',num2str(CCL21_block),'\n']);
             
            
             [p,x,cross_pairs,h] = cfd_FiberGenerator_flow_Sep26th(p,kk);
             if round(p.CollagenDensity*100)/100 ~= p.rho
                 p.n = floor(p.n/p.CollagenDensity*p.rho);
                 [p,x,cross_pairs,h] = cfd_FiberGenerator_flow_Sep26th(p,kk);
             end
             
             fprintf(['collagen density: ',num2str(p.CollagenDensity),'\n']);
             fprintf(['collagen fiber number: ',num2str(p.n),'\n']);
             x_cell = cell_initiation_flow_May17th(p,kk);
             i = 1;
             h = 0;
            
             F = TractionForceCalculator(p,x,x_cell);
             F = ProtursionForceCalculator(F,p,x,x_cell);
             MSD = zeros();
             movement = zeros();
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
            dlmwrite(['./',p.filename,'/',p.filename,'_movement.csv'], [], 'delimiter', ',');
            
            
            
                    
            x_cell.location_x_previous = 0;
            x_cell.location_y_previous = 0;
            x_cell.location_z_previous = 0;
            % hold on;
            for i = 1:p.tnind
                       %x_cell.v = zeros();
                for m = 1:100
                    for n = 1:3
                if abs(drag(m,n) + F.protrusion(m,n) + F.traction(m,n)) <= abs(F.frictionECM(m,n))
                   x_cell.v (m,n) = 0;
                else
                x_cell.v (m,n) = (drag (m,n) + F.protrusion (m,n) + F.traction (m,n) + F.frictionECM(m,n))./(p.viscosity_ECM*6*pi*7.5);
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
                movement = sqrt((x_cell.location_x - x_cell.location_x_source).^2 + (x_cell.location_y - x_cell.location_y_source).^2 + (x_cell.location_z - x_cell.location_z_source).^2);
                speed_magnitude = sqrt((x_cell.v(:,1)).^2 + (x_cell.v(:,2)).^2 + (x_cell.v(:,3)).^2);

                F = TractionForceCalculator(p,x,x_cell);
                F = ProtursionForceCalculator(F,p,x,x_cell);

                v_theta = acos(x_cell.v(:,3)*p.dt./sqrt((x_cell.v(:,1)*p.dt).^2 + (x_cell.v(:,2)*p.dt).^2 + (x_cell.v(:,3)*p.dt).^2))*180/pi;
                v_phi = atan2(x_cell.v(:,2)*p.dt,x_cell.v(:,1)*p.dt)*180/pi;

                tract_theta = acos(F.traction(:,3)./sqrt((F.traction(:,1)).^2 + (F.traction(:,2)).^2 + (F.traction(:,3)).^2))*180/pi;
                tract_phi = atan2(F.traction(:,2),F.traction(:,1))*180/pi;

                cell_theta = acos(x_cell.location_z./sqrt((x_cell.location_x).^2 + (x_cell.location_y).^2 + (x_cell.location_z).^2))*180/pi;
                cell_phi = atan2(x_cell.location_y,x_cell.location_x)*180/pi;
                
% %                 if mod(p.tnind,8) == 0
% %                     temp_theta = acos((x_cell.location_z - x_cell.location_z_previous)./sqrt((x_cell.location_x - x_cell.location_x_previous).^2 + (x_cell.location_y - x_cell.location_y_previous).^2 + (x_cell.location_z - x_cell.location_z_previous).^2))*180/pi;
% %                     streamline = (sum(temp_theta <= 45 | temp_theta >= 135) - sum(temp_theta > 45 & temp_theta < 135))/size(temp_theta,1);
% %                     directional = (sum(temp_theta <= 45) - sum(temp_theta >= 135))/(sum(temp_theta <= 45) + sum(temp_theta >= 135));
% %                     x_cell.location_x_previous = x_cell.location_x;
% %                     x_cell.location_y_previous = x_cell.location_y;
% %                     x_cell.location_z_previous = x_cell.location_z;
% %                     if CCL21_block == 0
% %                         if k == 0
% %                             dlmwrite('cfd_00um_streamline_control.csv', streamline,'delimiter',',','-append');
% %                             dlmwrite('cfd_00um_directional_control.csv', directional,'delimiter',',','-append');
% %                         elseif k == 1
% %                             dlmwrite('cfd_03um_streamline_control.csv', streamline,'delimiter',',','-append');
% %                             dlmwrite('cfd_03um_directional_control.csv', directional,'delimiter',',','-append');
% %                         elseif k == 2
% %                             dlmwrite('cfd_3um_streamline_control.csv', streamline,'delimiter',',','-append');
% %                             dlmwrite('cfd_3um_directional_control.csv', directional,'delimiter',',','-append');
% %                         end
% %                     elseif CCL21_block == 1
% %                         if k == 0
% %                             dlmwrite('cfd_00um_streamline_CCL21block.csv', streamline,'delimiter',',','-append');
% %                             dlmwrite('cfd_00um_directional_CCL21block.csv', directional,'delimiter',',','-append');
% %                         elseif k == 1
% %                             dlmwrite('cfd_03um_streamline_CCL21block.csv', streamline,'delimiter',',','-append');
% %                             dlmwrite('cfd_03um_directional_CCL21block.csv', directional,'delimiter',',','-append');
% %                         elseif k == 2
% %                             dlmwrite('cfd_3um_streamline_CCL21block.csv', streamline,'delimiter',',','-append');
% %                             dlmwrite('cfd_3um_directional_CCL21block.csv', directional,'delimiter',',','-append');
% %                         end                   
% %                     
% %                     end
% %                 end

                if i == p.tnind
                    temp_theta = acos((x_cell.location_z)./sqrt((x_cell.location_x).^2 + (x_cell.location_y).^2 + (x_cell.location_z).^2))*180/pi;
                    streamline = (sum(temp_theta <= 45 | temp_theta >= 135) - sum(temp_theta > 45 & temp_theta < 135))/size(temp_theta,1);
                    directional = (sum(temp_theta <= 45) - sum(temp_theta >= 135))/(sum(temp_theta <= 45) + sum(temp_theta >= 135));
                    
                    if CCL21_block == 0
                        if k == 0
                            dlmwrite('cfd_00um_streamline_control.csv', streamline,'delimiter',',','-append');
                            dlmwrite('cfd_00um_directional_control.csv', directional,'delimiter',',','-append');
                        elseif k == 1
                            dlmwrite('cfd_03um_streamline_control.csv', streamline,'delimiter',',','-append');
                            dlmwrite('cfd_03um_directional_control.csv', directional,'delimiter',',','-append');
                        elseif k == 2
                            dlmwrite('cfd_3um_streamline_control.csv', streamline,'delimiter',',','-append');
                            dlmwrite('cfd_3um_directional_control.csv', directional,'delimiter',',','-append');
                        end
                    elseif CCL21_block == 1
                        if k == 0
                            dlmwrite('cfd_00um_streamline_CCL21block.csv', streamline,'delimiter',',','-append');
                            dlmwrite('cfd_00um_directional_CCL21block.csv', directional,'delimiter',',','-append');
                        elseif k == 1
                            dlmwrite('cfd_03um_streamline_CCL21block.csv', streamline,'delimiter',',','-append');
                            dlmwrite('cfd_03um_directional_CCL21block.csv', directional,'delimiter',',','-append');
                        elseif k == 2
                            dlmwrite('cfd_3um_streamline_CCL21block.csv', streamline,'delimiter',',','-append');
                            dlmwrite('cfd_3um_directional_CCL21block.csv', directional,'delimiter',',','-append');
                        end                   
                    
                    end
                end

                
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
               % dlmwrite(['./',p.filename,'/',p.filename,'_R_fiberN.csv'], F.R_fiberN,'delimiter',',','-append');


                dlmwrite(['./',p.filename,'/',p.filename,'_Tx.csv'], F.traction(:,1)'*1e12, 'delimiter', ',','-append');
                dlmwrite(['./',p.filename,'/',p.filename,'_Ty.csv'], F.traction(:,2)'*1e12, 'delimiter', ',','-append');
                dlmwrite(['./',p.filename,'/',p.filename,'_Tz.csv'], F.traction(:,3)'*1e12, 'delimiter', ',','-append');
                dlmwrite(['./',p.filename,'/',p.filename,'_movement.csv'], movement', 'delimiter', ',','-append');

            end


            

            % speed output section
            speed = csvread(['./',p.filename,'/',p.filename,'_v.csv']);
            
             % last position direction output section
            [N,~] = histc(cell_theta,0:10:180);
            

            % MSD calculation !!
            cell_positions_x(p.tnind+1,:) = x_cell.location_x';
            cell_positions_y(p.tnind+1,:) = x_cell.location_y';
            cell_positions_z(p.tnind+1,:) = x_cell.location_z';
            MSD_complex = zeros(p.tnind,1);
            for i = 1:p.tnind
                MSD_complex(i) = sum(sum((cell_positions_x(i+1:end,:) - cell_positions_x(1:end-i,:)).^2 + (cell_positions_y(i+1:end,:) - cell_positions_y(1:end-i,:)).^2 +(cell_positions_z(i+1:end,:) - cell_positions_z(1:end-i,:)).^2))/(size(cell_positions_x(i+1:end,:),1)*size(cell_positions_x(i+1:end,:),2))/30/30;
            end
            dlmwrite(['./',p.filename,'/',p.filename,'_MSD_complex.csv'], MSD_complex, 'delimiter', ',');

            if j == 1
                figure('PaperPosition',[.25 .25 16 12]);
                loglog((1:p.tnind)*2,MSD_complex,'k-','linewidth',2);
                set(gca,'fontsize',12)
                title([num2str(round(p.CollagenDensity)),'mg/ml collagen density ',num2str(p.flow_v_z),'um/s interstitial flow'], 'Fontsize', 12);
                ylabel('MSD (um^2/hr^2)', 'Fontsize', 12);
                xlabel('time(min)', 'Fontsize', 12);
                set(gca,'linewidth',2)
                print('-dpng','-r400',['./figure/',p.filename,'_MSD_complex.png']);
            end

            logx = log10((1:p.tnind)*2)';
            logMSD = log10(MSD_complex);
            myfittype = fittype('a*x+b','independent','x');
            for i = 2:p.tnind
                [flogMSD,gof] = fit(logx(max([i - 20,1]):i),logMSD(max([i - 20,1]):i), myfittype);
                flogMSDrollingslope(i - 1) = flogMSD.a;
                flogMSDrollingintercept(i - 1) = flogMSD.b;
                if i <= p.tnind - 1 && i >= 10
                    [flogMSD,gof] = fit(logx(end - i:end),logMSD(end - i:end), myfittype);
                    flogMSDslope(p.tnind - i) = flogMSD.a;
                    flogMSDintercept(p.tnind - i) = flogMSD.b;
                end
            end

    %         overall_MSD_slope(l) = flogMSDrollingslope(100);
    %         overall_MSD_slope2(l) = flogMSDslope(100);


            mean_speed = mean(speed');
            std_speed = std(speed');

            figure('PaperPosition',[.25 .25 16 12]);
            errorbar((1:p.tnind)*2,mean_speed,std_speed);
            hold on;
            plot((1:p.tnind)*2,mean_speed,'r-','linewidth',2);
            set(gca,'fontsize',12)
            xlim([0 p.tnind*2]);
    %         ylim([0 0.2]);
            title([num2str(round(p.CollagenDensity)),'mg/ml collagen density ',num2str(p.flow_v_z),'um/s interstitial flow'], 'Fontsize', 12);
            ylabel('speed magnitude (um/h)', 'Fontsize', 12);
            xlabel('time(min)', 'Fontsize', 12);
            set(gca,'linewidth',2)
            print('-dpng','-r400',['./figure/',p.filename,'_speed_magnitude.png']);

            disp_z = csvread(['./',p.filename,'/',p.filename,'_z.csv']);
            figure('PaperPosition',[.25 .25 16 12]);
            errorbar((0:p.tnind)*2,mean(disp_z'),std(disp_z'));
            hold on;
            plot((0:p.tnind)*2,mean(disp_z'),'r-','linewidth',2);
            set(gca,'fontsize',12)
            xlim([0 p.tnind*2]);
            title([num2str(round(p.CollagenDensity)),'mg/ml collagen density ',num2str(p.flow_v_z),'um/s interstitial flow'], 'Fontsize', 12);
            ylabel('displacement (um)', 'Fontsize', 12);
            xlabel('time(min)', 'Fontsize', 12);
            set(gca,'linewidth',2)
            print('-dpng','-r400',['./figure/',p.filename,'_z_displacement.png']);




            l = l + 1;
         end
     end
     
 end
 
 
directional_control = csvread("cfd_03um_directional_control.csv");
streamline_control = csvread("cfd_03um_streamline_control.csv");
directional_CCL21block = csvread("cfd_03um_directional_CCL21block.csv");
streamline_CCL21block = csvread("cfd_03um_streamline_CCL21block.csv");
directional_control_mean = mean(directional_control);
streamline_control_mean = mean(streamline_control);
directional_CCL21block_mean = mean(directional_CCL21block);
streamline_CCL21block_mean = mean(streamline_CCL21block);
 
