
 clear;
 clc;
 close all;


l = 1;
kk = 2;
k = 0;
for jj = 1:4
    if jj == 1
        CSCG_size = 15;
    elseif jj == 2
        CSCG_size = 30;
    elseif jj == 3
        CSCG_size = 45;
    elseif jj == 4
        CSCG_size = 60;
    end
 for j = 1:1
     close all;
    
     filename = ['cdf_test_20um_',num2str(j)];
     clearvars p x cross_pairs;
     p = cfd_initiation_flow_May17th_02mg(filename,k,kk);
     p.bottom_scale_x = -20; % the coordinate of the fiber network cube starts, such as (200,200,200) to (220,220,220)
     p.bottom_scale_y = -20; % the coordinate of the fiber network cube starts, such as (200,200,200) to (220,220,220)
     p.bottom_scale_z = -CSCG_size/2; % the coordinate of the fiber network cube starts, such as (200,200,200) to (220,220,220)
     p.domain_scale_z = CSCG_size;
     p.domain_scale_y = 40;
     p.domain_scale_x = 40;
     p.V = p.domain_scale_x*p.domain_scale_y*p.domain_scale_z;

     [p,x,cross_pairs,h] = cfd_FiberGenerator_flow_Sep26th(p,kk);
     if round(p.CollagenDensity*100)/100 ~= 1
         fprintf('regenerate fibernetwork\n');
         p.n = floor(p.n/p.CollagenDensity*p.rho);
         [p,x,cross_pairs,h] = cfd_FiberGenerator_flow_Sep26th(p,kk);
     end
     x_overall = x;
     p.bottom_scale_x = 20; % the coordinate of the fiber network cube starts, such as (200,200,200) to (220,220,220)
     p.bottom_scale_y = 20; % the coordinate of the fiber network cube starts, such as (200,200,200) to (220,220,220)
     p.bottom_scale_z = CSCG_size/2; % the coordinate of the fiber network cube starts, such as (200,200,200) to (220,220,220)
     p.domain_scale_z = 250;
     p.domain_scale_y = 40;
     p.domain_scale_x = 40;
     
     [p,x,cross_pairs,h] = cfd_FiberGenerator_flow_Sep26th(p,1);
     if round(p.CollagenDensity*100)/100 ~= 2
         fprintf('regenerate fibernetwork\n');
         p.n = floor(p.n/p.CollagenDensity*p.rho);
         [p,x,cross_pairs,h] = cfd_FiberGenerator_flow_Sep26th(p,1);
     end
     
     
     
     x_cell = cell_initiation_flow_May17th(p,kk);
     i = 1;
     h = 0;
    %  h = cfd_output_flow_May18th2017(p,x,x_cell,1,h);
     fprintf(['alignment: ',num2str(kk),'\n']);
     fprintf(['flowspeed: ',num2str(k),'\n']);
     fprintf(['trial: ',num2str(j),'\n']);
    %  x.rotcenter = p.rotcenter;
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

    %     for j = 1:x_cell.N
    %         %plot3(x_cell.location_x,x_cell.location_y,x_cell.location_z,'r.','Markersize',12);
    %         hold on;
    %         %c = [(j - 1)*1/(x_cell.N-1) 1-(j - 1)*1/(x_cell.N-1) (j - 1)*1/(x_cell.N-1)];
    %         c = [1-(j - 1)*1/(x_cell.N-1) (j - 1)*1/(x_cell.N-1) (j - 1)*1/(x_cell.N-1)];
    %         quiver3(original_location_x(j),original_location_y(j),original_location_z(j),x_cell.v(j,1)*p.dt,x_cell.v(j,2)*p.dt,x_cell.v(j,3)*p.dt,'color',c,'AutoScale','off','ShowArrowHead','off');
    %         hold on;
    %     end
        %pause();
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


        dlmwrite(['./',p.filename,'/',p.filename,'_Tx.csv'], F.traction(:,1)', 'delimiter', ',','-append');
        dlmwrite(['./',p.filename,'/',p.filename,'_Ty.csv'], F.traction(:,2)', 'delimiter', ',','-append');
        dlmwrite(['./',p.filename,'/',p.filename,'_Tz.csv'], F.traction(:,3)', 'delimiter', ',','-append');
        dlmwrite(['./',p.filename,'/',p.filename,'_movement.csv'], movement', 'delimiter', ',','-append');

    end


    % last position output section       
    overall_last_location_z_mean(l) = mean(x_cell.location_z);
    overall_last_location_z_std(l) = std(x_cell.location_z);


    overall_last_location_x_mean(l) = mean(x_cell.location_x);
    overall_last_location_x_std(l) = std(x_cell.location_x);


    overall_last_location_y_mean(l) = mean(x_cell.location_y);
    overall_last_location_y_std(l) = std(x_cell.location_y);



    % speed output section
    speed = csvread(['./',p.filename,'/',p.filename,'_v.csv']);
    overall_speed_mean(l) = mean(speed(:));
    overall_speed_std(l) = std(speed(:));

     % last position direction output section
    [N,~] = histc(cell_theta,0:10:180);
    overall_theta(l,:) = N;

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
        loglog((1:p.tnind)*p.dt/60,MSD_complex,'k-','linewidth',2);
        set(gca,'fontsize',12)
        title([num2str(round(p.CollagenDensity)),'mg/ml collagen density ',num2str(p.flow_v_z),'um/s interstitial flow'], 'Fontsize', 12);
        ylabel('MSD (um^2/hr^2)', 'Fontsize', 12);
        xlabel('time(min)', 'Fontsize', 12);
        set(gca,'linewidth',2)
        print('-dpng','-r400',['./figure/',p.filename,'_MSD_complex.png']);
    end

    logx = log10((1:p.tnind)*p.dt/60)';
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
%         

    mean_speed = mean(speed');
    std_speed = std(speed');

    figure('PaperPosition',[.25 .25 16 12]);
    errorbar((1:p.tnind)*p.dt/60,mean_speed,std_speed);
    hold on;
    plot((1:p.tnind)*p.dt/60,mean_speed,'r-','linewidth',2);
    set(gca,'fontsize',12)
    xlim([0 p.tnind*p.dt/60]);
%         ylim([0 0.2]);
    title([num2str(round(p.CollagenDensity)),'mg/ml collagen density ',num2str(p.flow_v_z),'um/s interstitial flow'], 'Fontsize', 12);
    ylabel('speed magnitude (um/h)', 'Fontsize', 12);
    xlabel('time(min)', 'Fontsize', 12);
    set(gca,'linewidth',2)
    print('-dpng','-r400',['./figure/',p.filename,'_speed_magnitude.png']);

    disp_z = csvread(['./',p.filename,'/',p.filename,'_z.csv']);
    figure('PaperPosition',[.25 .25 16 12]);
    errorbar((0:p.tnind)*p.dt/60,mean(disp_z'),std(disp_z'));
    hold on;
    plot((0:p.tnind)*p.dt/60,mean(disp_z'),'r-','linewidth',2);
    set(gca,'fontsize',12)
    xlim([0 p.tnind*p.dt/60]);
    title([num2str(round(p.CollagenDensity)),'mg/ml collagen density ',num2str(p.flow_v_z),'um/s interstitial flow'], 'Fontsize', 12);
    ylabel('displacement (um)', 'Fontsize', 12);
    xlabel('time(min)', 'Fontsize', 12);
    set(gca,'linewidth',2)
    print('-dpng','-r400',['./figure/',p.filename,'_z_displacement.png']);

 end 
end
        