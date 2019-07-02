function F = ProtursionForceCalculator(F,p,x,x_cell)


index_tmp = rand(x_cell.N,1);
F.protursion_theta = p.temp_pool(floor((length(p.temp_pool) - 1).*index_tmp)+1);
F.protursion_phi = 360*rand(x_cell.N,1);


F.protrusion = [p.mag_protursion*sin(F.protursion_theta/180*pi).*cos(F.protursion_phi/180*pi) p.mag_protursion*sin(F.protursion_theta/180*pi).*sin(F.protursion_phi/180*pi) p.mag_protursion*cos(F.protursion_theta/180*pi)];
