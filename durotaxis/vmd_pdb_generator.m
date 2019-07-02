% This is the file to convert all one-time cordinates to filename.pdb
% input      
%           tind --- the time index, or the frame index that starts with 1
%  x_seg_cor_glo --- the x cordinates of each node                 
%  y_seg_cor_glo --- the y cordinates of each node
%  z_seg_cor_glo --- the z cordinates of each node
%       filenmae --- the filename.xyz as the output
% output  
%        a .pdb file in text format as the input of the vmd containing the
%        frame tind all nodes coordinates.
% Meng Sun
% Muhammad Zaman Lab
% Biomedical Engineering
% Boston University
% email: ttsunmeng@gmail.com
                     
function vmd_pdb_generator(tind,fiber_nodes_num,x_seg_cor_glo,y_seg_cor_glo,z_seg_cor_glo,filename)



N = sum(fiber_nodes_num);%


if tind == 1
    fileID = fopen(['./',filename,'/',filename,'.pdb'],'w');
    fprintf(fileID,'%6s\n','REMARK');
else
    fileID = fopen(['./',filename,'/',filename,'.pdb'],'a');
end
    fprintf(fileID,'%5s %8d\n','MODEL',tind);
for i = 1:N
    fprintf(fileID,'%6s %4d %17s %6.3f %5.3f %5.3f %18s\n','ATOM  ',i,' ',x_seg_cor_glo(i),y_seg_cor_glo(i),z_seg_cor_glo(i),'C');
end
fprintf(fileID,'%6s\n','ENDMDL');
fclose(fileID);