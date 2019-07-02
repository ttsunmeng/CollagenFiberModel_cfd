function vmd_psf_generator(fiber_nodes_num,filename)
N = sum(fiber_nodes_num);%
fileID = fopen(['./',filename,'/',filename,'.psf'],'w');
fprintf(fileID,'%8s\n','PSF CMAP');
fprintf(fileID,'\n');
fprintf(fileID,'%8d %7s\n',7,'!NTITLE');
fprintf(fileID,'\n');
fprintf(fileID,'%8d %6s\n',N,'!NATOM');
for i = 1:N
    fprintf(fileID,'%8d %1s %4d %53s\n',i,'C',1,'   LYS  C      11    0.00000        1.0000          0');
end

n = length(fiber_nodes_num);%
% a = zeros(N,1);%
b(2:n + 1) = cumsum(fiber_nodes_num);%b is the previous node of each left node.
b(1) = 0;
c = zeros(N - n,1);
d = zeros(N - n,1);
for i = 1:n
%     a(b(i) + 1:b(i + 1)) = 1:fiber_nodes_num(i);
    c(b(i) + 2 - i:b(i + 1) - i) = b(i) + 1:b(i) + fiber_nodes_num(i) - 1;
    d(b(i) + 2 - i:b(i + 1) - i) = b(i) + 2:b(i) + fiber_nodes_num(i);
end

bond_xyz = [c';d'];
fprintf(fileID,'\n');
fprintf(fileID,'%8d %13s\n',length(bond_xyz),'!NBOND: bonds');
for i = 1:4:fix(length(bond_xyz(1,:))/4)*4
    fprintf(fileID,' %7d %7d %7d %7d %7d %7d %7d %7d\n',transpose(bond_xyz(:,i)),transpose(bond_xyz(:,i+1)),transpose(bond_xyz(:,i+2)),transpose(bond_xyz(:,i+3)));
end
if mod(length(bond_xyz(1,:)),4) == 1
    fprintf(fileID,' %7d %7d \n',transpose(bond_xyz(:,i)),transpose(bond_xyz(:,i+1)));
elseif mod(length(bond_xyz(1,:)),4) == 2
    fprintf(fileID,' %7d %7d %7d %7d\n',transpose(bond_xyz(:,i)),transpose(bond_xyz(:,i+1)),transpose(bond_xyz(:,i+2)));
elseif mod(length(bond_xyz(1,:)),4) == 3
    fprintf(fileID,' %7d %7d %7d %7d %7d %7d\n',transpose(bond_xyz(:,i)),transpose(bond_xyz(:,i+1)),transpose(bond_xyz(:,i+2)),transpose(bond_xyz(:,i+3)));
end

fclose(fileID);
