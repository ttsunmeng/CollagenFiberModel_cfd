function x_cell = cell_initiation_flow_May17th(p,kk)

x_cell.N = 100; %cell number

x_cell.R = ones(x_cell.N,1)*7.5; % cell radius in m.


x_cell.location_x = p.bottom_scale_x + p.domain_scale_x*rand(x_cell.N,1);
x_cell.location_y = p.bottom_scale_y + p.domain_scale_y*rand(x_cell.N,1);
x_cell.location_z = p.bottom_scale_z + p.domain_scale_z*rand(x_cell.N,1);

x_cell.location_x_source = x_cell.location_x;
x_cell.location_y_source = x_cell.location_y;
x_cell.location_z_source = x_cell.location_z;
x_cell.v = zeros(x_cell.N,3);
x_cell.cubindex_x = 0;
x_cell.cubindex_y = 0;
x_cell.cubindex_z = 0;