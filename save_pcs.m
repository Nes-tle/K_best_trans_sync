function [] = save_pcs(pc, Results, foldername)
%
for sId = 1 : size(pc, 1)
    surf = reshape(pc(sId,:,:), [8192,3])';
    R_opt = Results{sId}.R_opts{1};
    t_opt = Results{sId}.t_opts{1};
    surf = R_opt'*surf - R_opt'*t_opt*ones(1, size(surf,2));
    str = sprintf('%d.xyz', sId);
    f_id = fopen([foldername, str], 'w');
    for i = 1 : size(surf,2)
        fprintf(f_id, '%f %f %f\n', surf(1, i), surf(2, i), surf(3,i));
    end
    fclose(f_id);
end