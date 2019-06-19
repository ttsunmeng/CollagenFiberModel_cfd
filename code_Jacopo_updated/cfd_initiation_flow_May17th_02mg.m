function p = cfd_initiation_flow_May17th_02mg(filename,k,kk)
    p.filename = filename;
    % The initiation of the basic domains
    p.bottom_scale = -50; % the coordinate of the fiber network cube starts, such as (200,200,200) to (220,220,220)
    
    p.domain_scale_z = 100;
    p.domain_scale_y = 100;
    p.domain_scale_x = 100;
    p.n = 3.3118e+04;
    
    p.V = p.domain_scale_x*p.domain_scale_y*p.domain_scale_z;
    
    if k ~= 0
        if kk == 1
        flow_list = [0.3,0.5,1,1.5,2,2.5,3];
        p.flow_v = [0 0 flow_list(k)]; % um/s from Polacheck et al 2011 PNAS.
        p.flow_v_z = flow_list(k);
        chemotaxis_table = [0.586,0.623,0.657,0.677,0.691,0.7,0.706;0.431,0.412,0.39,0.37,0.356,0.355,0.333;0.218,0.183,0.152,0.132,0.12,0.115,0.103;0.166,0.133,0.106,0.091,0.082,0.079,0.071;0.142,0.114,0.092,0.079,0.071,0.068,0.062];
        flow_v_table = [0.3;0.5;1;1.5;2;2.5;3];
        theta_table = [0;45;90;135;180];
        [X,Y] = meshgrid(flow_v_table,theta_table);
        [Xq,Yq] = meshgrid(norm(p.flow_v),(0:1:180)');
        Vq = 100*interp2(X,Y,chemotaxis_table,Xq,Yq);
        temp_pool = [];
        for i = 1:181
            temp_pool = [temp_pool;ones(floor(Vq(i)),1)*(i-1)];
        end
        p.temp_pool = temp_pool;
        elseif kk == 2
        flow_list = [0.3,0.5,1,1.5,2,2.5,3];
        p.flow_v = [0 0 flow_list(k)]; % um/s from Polacheck et al 2011 PNAS.
        p.flow_v_z = flow_list(k);
        chemotaxis_table = [0.53,0.559,0.594,0.615,0.63,0.64,0.654;0.415,0.389,0.357,0.33,0.308,0.291,0.28;0.243,0.208,0.176,0.155,0.14,0.129,0.121;0.205,0.172,0.142,0.124,0.112,0.103,0.097;0.18,0.151,0.125,0.11,0.10,0.092,0.087];
        flow_v_table = [0.3;0.5;1;1.5;2;2.5;3];
        theta_table = [0;45;90;135;180];
        [X,Y] = meshgrid(flow_v_table,theta_table);
        [Xq,Yq] = meshgrid(norm(p.flow_v),(0:1:180)');
        Vq = 100*interp2(X,Y,chemotaxis_table,Xq,Yq);
        temp_pool = [];
        for i = 1:181
            temp_pool = [temp_pool;ones(floor(Vq(i)),1)*(i-1)];
        end
        p.temp_pool = temp_pool;
            
        elseif kk == 3   
        flow_list = [0.3,0.5,1,1.5,2,2.5,3];
        p.flow_v = [0 0 flow_list(k)]; % um/s from Polacheck et al 2011 PNAS.
        p.flow_v_z = flow_list(k);
        chemotaxis_table = [0.51,0.533,0.566,0.588,0.603,0.62,0.627;0.433,0.414,0.391,0.367,0.349,0.332,0.316;0.242,0.204,0.164,0.14,0.124,0.112,0.102;0.211,0.18,0.148,0.13,0.118,0.109,0.102;0.211,0.187,0.162,0.147,0.136,0.127,0.12];
        flow_v_table = [0.3;0.5;1;1.5;2;2.5;3];
        theta_table = [0;45;90;135;180];
        [X,Y] = meshgrid(flow_v_table,theta_table);
        [Xq,Yq] = meshgrid(norm(p.flow_v),(0:1:180)');
        Vq = 100*interp2(X,Y,chemotaxis_table,Xq,Yq);
        temp_pool = [];
        for i = 1:181
            temp_pool = [temp_pool;ones(floor(Vq(i)),1)*(i-1)];
        end
        p.temp_pool = temp_pool;       
        end
    elseif k == 0
        p.flow_v = [0 0 0];
        p.flow_v_z = 0;
        p.temp_pool = (0:180)';
    end
    
    
    p.mag_protursion = 50e-12; % 50-500pN: half goes to protursion (50pN) and half to traction?
    p.E = 1e-4; % 100kPa = 10^(-7)N/um^2. 100MPa = 1e-4N/um^2
    p.ctraction = 3e-8; % So that the traction force on each fiber is 3pN = p.ctraction*p.E.1-100pN per bond
    
    p.cresistance = 177*3e-4;%15/3600*p.dt*p.n/p.V; % So the resistance force we assume is 10 times smaller here.
    % Effective area of the cell is pi*R*R=177um^2, U = 15um/h =
    % 15/3600*p.dt
    
    
    p.rho = 2; % collagen concentration 2mg/ml;
    %p.n = floor(2500/3.3/8000*p.rho*p.V);
    
    %the number of collagen fibers/fibrils
    p.fiber_length_a = 0.971750840346117;                                      % from 4mg/ml FIRE length data estimated by gamma distribution.
    p.fiber_length_b = 1.748381900546554;                                      % from 4mg/ml FIRE length data estimated by gamma distribution.
    p.fiber_length_base = 1.9859;
    p.fiber_length_mean = p.fiber_length_a*p.fiber_length_b + p.fiber_length_base;
    p.persistent_length = 1.0;
    
    
    AI = 0.9;
    checkTable = [0,1;5,0.987;10,0.957;15,0.91;20,0.843;25,0.77;30,0.632;35,0.52;40,0.43;45,0.301;50,0.216;55,0.133;60,0.12;65,0.0958;70,0.0789;180,0];
    p.theta_std = interp1(checkTable(:,2),checkTable(:,1),AI);
    % These data are from Sivakumar et al, 2010, The influence of discoidin 
    % domain receptor 2 on the persistence length of collagen type I fibers
    p.seg_diameter = 0.42; %ranges a lot, but can choose the average of the scale to be 100nm from Sivakumar et al, 2010, The influence of discoidin domain receptor 2 on the persistence length of collagen type I fibers
    % Stein et al 2008 (DOI: 10.1111/j.1365-2818.2008.02141.x) also used 30nm. 
    p.A = pi*p.seg_diameter*p.seg_diameter/4;
    p.Irot = 1/64*pi*p.seg_diameter*p.seg_diameter*p.seg_diameter*p.seg_diameter; % Second moment of inertia for circular cross section for all cartesian axes
    p.kbend = p.E*p.Irot/p.fiber_length_mean;
    %     p.orientation_random = 1; % if 1, it is all random.
%     p.orientation_theta = 0;
%     p.orientation_theta_std = 0;
%     p.orientation_phi = 0;
%     p.orientation_phi_std = 0;

    p.massdensity = 1300/10^12; % unit in 1.3 g/cm3 in 2013 Fibrillar structure and elasticity of hydrating collagen- A quantitative multiscale approach.
    % 1.3 g/cm3 = 1300mg/ml;
    
    
    p.viscosity = 1e-15; %1e-9/4; %unit in 0.001 Pa*s = N*s/m^2 = 0.001*10^(-12) N*s/um^2 = 10^(-15) This is the water. 37��C water is 0.0007 Pa*s 
    % Here we use the zero shear viscosity rather than the interstitial flow
    % viscosity to make the force magnitude resonable: 10^5 Pa*s = 10^(-7)
    % N*s/um^2 This is overall collagen gel
   
    % fomula from Kim et al, 2009, Computational Analysis of Viscoelastic Properties of Crosslinked Actin Networks
    p.viscosity_ECM = 1e-10; %1e-9/4;
    
% %     %crosslink data
% %     p.distance_crosslinking = 0.5; 
% %     % Head et al, 2003 gives that spacing between apparent crosslinks is
% %     % Possion distribution, with <Lc> = 2.5um ?Stein et al, 2008?.
% %     p.break_force = 2.34e-1; % assume 10^(-6) N;
% %     % To break one fibril into monomers is about 10^(-9) N = 1.5nN.
% %     % mechanics properties of fibrils initiation%%%%%%%%%%%%%%%%%%
% %     p.possible_crosslinks = 1;
% %     
    
    %% time evolution initiation
    p.Solver_option = 2; % if 1, the Euler Solver; if 2, the Runge Kuta Solver;
    if p.Solver_option == 1
        p.dt = 1e-5; 
       
    elseif p.Solver_option == 2
        p.dt = 120; % time interval 2 min
        
    end    
    
    
    
    
    p.Tend = 75*3600;  
    p.tnind = p.Tend/p.dt;
    

    
    %% the setup of the display    
    p.vis_in_matlab = 1; % Whether to visualize the simulation in matlab directly.
    p.safe_vis = 1;
    p.vis_in_VMD = 0; % Whether to generate the vmd files.
    p.output_orientation = 1;
    p.output_selected = 1; % Whether to generate selected node positions only in csv files with the stress.
    p.output_all = 0; % Whether to generate each node position in csv files.
    p.output_precalculated_strain_in_excel_middle = 1;
    
    mkdir(p.filename);
    rmdir(p.filename,'s');
    mkdir(p.filename);

    
