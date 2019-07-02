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
        flow_list = [0.3,0.5,1,1.5,2,2.5,3];
        p.flow_v = [0 0 flow_list(k)]; % um/s from Polacheck et al 2011 PNAS.
        p.flow_v_z = flow_list(k);
        chemotaxis_table = [0.55,0.56,0.57,0.59,0.59,0.6,0.6;0.435,0.45,0.36,0.305,0.28,0.265,0.255;0.31,0.285,0.225,0.185,0.17,0.155,0.145;0.26,0.23,0.185,0.15,0.135,0.125,0.115;0.23,0.2,0.16,0.13,0.12,0.11,0.1];
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
    elseif k == 0
        p.flow_v = [0 0 0];
        p.flow_v_z = 0;
        p.temp_pool = (0:180)';
    end
    
    
    p.mag_protursion = 5e-11; % 50-500pN: half goes to protursion (50pN) and half to traction?
    p.E = 1e-7; % 100kPa = 10^(-7)N/um^2.
    p.ctraction = 10e-6; % So that the traction force on each fiber is 3pN = p.ctraction*p.E.1-100pN per bond
    
    
    p.cresistance = 1.77e2*300; % So the resistance force we assume is 10 times smaller here.Effective area of the cell is pi*R*R=177um^2
    
    p.rho = 2; % collagen concentration 2mg/ml;
    %p.n = floor(2500/3.3/8000*p.rho*p.V);
    
    %the number of collagen fibers/fibrils
    p.fiber_length_mean = 14; % length parameters micrometers
    p.fiber_length_std = 7;
    p.fiber_length_lowerlimit = 5;
    p.fiber_length_upperlimit = 23;
    
    
    AI = 0.9;
    checkTable = [0,1;5,0.987;10,0.957;15,0.91;20,0.843;25,0.77;30,0.632;35,0.52;40,0.43;45,0.301;50,0.216;55,0.133;60,0.12;65,0.0958;70,0.0789;180,0];
    p.theta_std = interp1(checkTable(:,2),checkTable(:,1),AI);
    % These data are from Sivakumar et al, 2010, The influence of discoidin 
    % domain receptor 2 on the persistence length of collagen type I fibers
    p.seg_diameter = 0.1; %ranges a lot, but can choose the average of the scale to be 100nm from Sivakumar et al, 2010, The influence of discoidin domain receptor 2 on the persistence length of collagen type I fibers
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
    
       
    p.tnind = 300;
    p.Tend = p.tnind*p.dt;  

    

    
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

    
