function [p,x,cross_pairs,h] = cfd_FiberGenerator_flow_Sep26th(p,kk)

    %% The lengths of fibers
    fiber_length = zeros(p.n,1);
    x_tmp = zeros(p.n,2);
    y = zeros(p.n,2);
    z = zeros(p.n,2);
    theta = zeros(p.n,1);
    phi = zeros(p.n,1);
    p.I = zeros(p.n,1);
    p.rotcenter = zeros(p.n,1);
    k = 1;
    if kk == 1
        theta = 90 + p.theta_std*randn(p.n,1); %p.theta_std*randn(p.n,1) %180*rand(p.n,1);
    elseif kk == 2
        theta = 90 + p.theta_std*randn(p.n,1); %p.theta_std*randn(p.n,1) %180*rand(p.n,1);
    elseif kk == 3
        theta = 180*rand(p.n,1); %p.theta_std*randn(p.n,1) %180*rand(p.n,1);
    elseif kk == 4
        theta = p.theta_std*randn(p.n,1); %p.theta_std*randn(p.n,1) %180*rand(p.n,1);
    end
    
    theta(theta<0) = -theta(theta<0);
    theta(theta > 360) = theta(theta > 360) - 360;
    theta(theta > 720) = theta(theta > 720) - 720;
    theta(theta > 180) = 360 - theta(theta > 180);
    theta(theta==180) = 0;
    while k <= p.n
        l = 1;
        while l == 1 || z(k,2) > p.domain_scale_z + p.bottom_scale_z || x_tmp(k,2) > p.domain_scale_x + p.bottom_scale_x || y(k,2) > p.domain_scale_y + p.bottom_scale_y...
                || z(k,1) < p.bottom_scale_z || x_tmp(k,1) < p.bottom_scale_x || y(k,1) < p.bottom_scale_y
            
%                 fiber_length(k) = gamrnd(p.fiber_length_a,p.fiber_length_b,1,1) + p.fiber_length_base;

                
                fiber_length(k) = p.fiber_length_mean + p.fiber_length_std*randn;


            if kk == 1
                phi(k) = p.theta_std*randn(1);
            elseif kk == 2
                phi(k) = 360*rand(1); % p.theta_std*randn(1);
            elseif kk == 3
                phi(k) = 360*rand(1); % p.theta_std*randn(1);
            elseif kk == 4
                phi(k) = 360*rand(1); % p.theta_std*randn(1);
            end


            x_tmp(k,1) = p.domain_scale_x*rand(1) + p.bottom_scale_x;
            y(k,1) = p.domain_scale_y*rand(1) + p.bottom_scale_y;
            z(k,1) = p.domain_scale_z*rand(1) + p.bottom_scale_z;
            
            z(k,2) = z(k,1) + fiber_length(k).*cos(theta(k)/180*pi);
            y(k,2) = y(k,1) + fiber_length(k).*sin(theta(k)/180*pi).*sin(phi(k)/180*pi);
            x_tmp(k,2) = x_tmp(k,1) + fiber_length(k).*sin(theta(k)/180*pi).*cos(phi(k)/180*pi);

            tmp = z(k,2);
            z(k,2) = z(k,1);
            z(k,1) = tmp;
            tmp = y(k,2);
            y(k,2) = y(k,1);
            y(k,1) = tmp;
            tmp = x_tmp(k,2);
            x_tmp(k,2) = x_tmp(k,1);
            x_tmp(k,1) = tmp;

            p.m(k) = p.massdensity*p.A*fiber_length(k);
            p.rotcenter(k) = rand;
            p.I(k) = (abs(0.5 - p.rotcenter(k))*fiber_length(k))^2*p.m(k) + 1/12*p.m(k)*fiber_length(k)^2; % units:[m]*um^2
            l = l + 1;
        end
        k = k + 1;
    end
    n = p.n;

    %% The tracking nodes of fibers and their coordinates.
    p.fiber_length = fiber_length;
    x.cor_glo = [x_tmp(:,1) + p.rotcenter.*(x_tmp(:,2) - x_tmp(:,1)) y(:,1) + p.rotcenter.*(y(:,2)-y(:,1)) z(:,1) + p.rotcenter.*(z(:,2)-z(:,1))];
    x.x_seg_start = x_tmp(:,1);
    x.x_seg_end = x_tmp(:,2);
    
    x.norm_x = sin(theta/180*pi).*cos(phi/180*pi);
    x.norm_y = sin(theta/180*pi).*sin(phi/180*pi);
    x.norm_z = cos(theta/180*pi);
    
    x.z_seg_start = z(:,1);
    x.z_seg_end = z(:,2);
    
    x.y_seg_start = y(:,1);
    x.y_seg_end = y(:,2);
    
    
%     x.start = [x.x_seg_start x.y_seg_start x.z_seg_start];
%     x.end = [x.x_seg_end x.y_seg_end x.z_seg_end];
    
    x.phi = phi;
    x.theta = theta;
    x.theta_new = x.theta;
    
    
    
%     tic
%     cross_pairs = crosslinkMarch3rd(p.distance_crosslinking, fiber_nodes_num,x.x_seg_cor_glo,x.y_seg_cor_glo,x.z_seg_cor_glo);      
%     toc
%     num_crx = size(cross_pairs,1);
%     display('the number of crosslinks and the ratio per fiber')
%     display([num_crx num_crx/p.n]);
    cross_pairs = [];
%     fprintf(['collagen density: ',num2str(sum(fiber_length)*p.A/(p.V+eps)*p.massdensity*10^12),'\n']);
%     fprintf(['collagen fiber number: ',num2str(p.n),'\n']);
    % pause();
    p.CollagenDensity = sum(fiber_length)*p.A/(p.V+eps)*p.massdensity*10^12;
    %% Output files
    h = 0;
    i = 1;
    
    p.n = n;
    