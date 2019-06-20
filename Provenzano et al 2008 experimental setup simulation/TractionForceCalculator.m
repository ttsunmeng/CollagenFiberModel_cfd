function F = TractionForceCalculator(p,x,x_cell)

if abs(x_cell.location_x) > p.general_scale/2
    x_cell.location_x = p.general_scale/2+sign(x_cell.location_x)*mod(abs(x_cell.location_x),p.general_scale/2);
    disp('out of scale in x');
end
if abs(x_cell.location_y) > p.general_scale/2
    x_cell.location_y = p.general_scale/2+sign(x_cell.location_y)*mod(abs(x_cell.location_y),p.general_scale/2);
    disp('out of scale in y');
end    
if abs(x_cell.location_z) > p.general_scale/2 + p.extension_scale
    error('out of scale in z');
end      

for i = 1:x_cell.N
    
    index_start = find(((x.x_seg_start - x_cell.location_x(i)).^2 + (x.y_seg_start - x_cell.location_y(i)).^2 + (x.z_seg_start - x_cell.location_z(i)).^2 < x_cell.R(i)^2)==1);
    index_end = find(((x.x_seg_end - x_cell.location_x(i)).^2 + (x.y_seg_end - x_cell.location_y(i)).^2 + (x.z_seg_end - x_cell.location_z(i)).^2 < x_cell.R(i)^2)==1);

    common_index = index_start(ismember(index_start,index_end));
    index_start_new = index_start;
    index_end_new = index_end;
    index_start_new(ismember(index_start,common_index)) = [];
    index_end_new(ismember(index_end,common_index)) = [];
    
    F.T_fiberN(i) = size(index_start_new,1) + size(index_end_new,1) + size(common_index,1);
    
   
    
    F.traction(i,:) = p.ctraction*p.E*[(sum(x.norm_x(index_start_new)) - sum(x.norm_x(index_end_new))) (sum(x.norm_y(index_start_new)) - sum(x.norm_y(index_end_new))) (sum(x.norm_z(index_start_new)) - sum(x.norm_z(index_end_new)))];
    
    if isempty(common_index)
        F.frictionECM(i,1) = 0;
        F.frictionECM(i,2) = 0;
        F.frictionECM(i,3) = 0;
        continue;
    end
    
    all_contact_index = [];
    project_vector = [];  
    
    project_vector = ([x.x_seg_start - x_cell.location_x(i),x.y_seg_start - x_cell.location_y(i),x.z_seg_start - x_cell.location_z(i)] - ...
        repmat(dot([x.x_seg_start - x_cell.location_x(i),x.y_seg_start - x_cell.location_y(i),x.z_seg_start - x_cell.location_z(i)],[x.x_seg_start - x.x_seg_end,x.y_seg_start - x.y_seg_end,x.z_seg_start - x.z_seg_end],2),1,3).*...
        [x.x_seg_start - x.x_seg_end,x.y_seg_start - x.y_seg_end,x.z_seg_start - x.z_seg_end]);
    
    % refer to
    % https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    % about vector formulation
    
    
    project_point = project_vector + repmat([x_cell.location_x(i),x_cell.location_y(i),x_cell.location_z(i)],p.n,1);

    
    all_contact_index = find(sqrt((x.x_seg_start - x_cell.location_x(i)).^2 + (x.y_seg_start - x_cell.location_y(i)).^2 + (x.z_seg_start - x_cell.location_z(i)).^2) <= x_cell.R(i) | ...
        sqrt((x.x_seg_end - x_cell.location_x(i)).^2 + (x.y_seg_end - x_cell.location_y(i)).^2 + (x.z_seg_end - x_cell.location_z(i)).^2) <= x_cell.R(i) |...
         (norm(project_vector) <= x_cell.R(i) & (dot(project_point - [x.x_seg_start,x.y_seg_start,x.z_seg_start], project_point - [x.x_seg_end,x.y_seg_end,x.z_seg_end],2) <= 0)));
   
    
    % the first term is the distance from cell center to fiber start point
    % the second term is the distance from cell center to fiber end point
    % the third term is the distance from cell center to the fiber line is
    % smaller and at the same time this projection point on this fiber.
    F.R_fiberN(i) = length(all_contact_index);
    dir_f_norm = [x.norm_x(all_contact_index) x.norm_y(all_contact_index) x.norm_z(all_contact_index)];
    x_cell_v_norm = x_cell.v(i,:)/(eps + realsqrt(sum(x_cell.v(i,:).^2,2)));
    dir_perpen_f_norm = cross(dir_f_norm,cross(dir_f_norm,repmat(x_cell_v_norm,length(all_contact_index),1)));
    
    cos_theta = sum(repmat(x_cell_v_norm,length(all_contact_index),1).*dir_f_norm,2);
    theta_bend = acos(abs(cos_theta));
    sin_theta = sqrt(1 - (cos_theta).^2);
    F.frictionECM(i,1) = p.cresistance*p.kbend*realsqrt(sum(x_cell.v(i,:).^2,2))*sum(theta_bend.*sin_theta.*dir_perpen_f_norm(:,1)./p.fiber_length(all_contact_index));
    F.frictionECM(i,2) = p.cresistance*p.kbend*realsqrt(sum(x_cell.v(i,:).^2,2))*sum(theta_bend.*sin_theta.*dir_perpen_f_norm(:,2)./p.fiber_length(all_contact_index));
    F.frictionECM(i,3) = p.cresistance*p.kbend*realsqrt(sum(x_cell.v(i,:).^2,2))*sum(theta_bend.*sin_theta.*dir_perpen_f_norm(:,3)./p.fiber_length(all_contact_index));
    %F.frictionECM(i,1) = p.cresistance*p.kbend*realsqrt(sum(x_cell.v(i,:).^2,2))*sum(sin_theta.*dir_perpen_f_norm(:,1)./14);
    %F.frictionECM(i,2) = p.cresistance*p.kbend*realsqrt(sum(x_cell.v(i,:).^2,2))*sum(sin_theta.*dir_perpen_f_norm(:,2)./14);
    %F.frictionECM(i,3) = p.cresistance*p.kbend*realsqrt(sum(x_cell.v(i,:).^2,2))*sum(sin_theta.*dir_perpen_f_norm(:,3)./14);


end