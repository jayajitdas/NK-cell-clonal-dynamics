function subsets = collect_subset_data(Ly49H_pos)

% Generates a output matrix with subset data.  Format of output is as
% follows:
% Rows are mice, columns are subset, 3rd dimension is time
% Subsets are [tom+27+ tom+DP tom+11b+ tom+DN tom-27+ tom-DP tom-11b+
% tom-DN]

days = {'1','3','7','14','23','29','35'}; %No D39 because some data doesn't exist
opts = detectImportOptions('D1_Bleed_Darrens_Gates.xls');
opts = setvartype(opts,'x_Cells','double');
for i=1:length(days)
    if Ly49H_pos
        input = readtable(strcat('D',days{i},'_Bleed_Darrens_Gates.xls'),opts);
        tot_cells = input.x_Cells(strcmp(input.Depth,'> > >'));
        my_subsets = input.x_Cells(strcmp(input.Depth,'> > > > > > >'));
    else
        input = readtable(strcat('./Ly49H-/D',days{i},'_Bleed_Darrens_Gates.xls'),opts);
        tot_cells = input.x_Cells(strcmp(input.Depth,'> > >'));
        my_subsets = input.x_Cells(contains(input.Name,'Ly49H-/'));
    end
    %check = input.Name(strcmp(input.Depth,'> > > > > > >'))
    clear tom_pos tom_neg
    for j=1:length(my_subsets)/8
        tom_pos(j,:) = my_subsets(j*8-7:j*8-4)/sum(my_subsets(j*8-7:j*8));
        tom_neg(j,:) = my_subsets(j*8-3:j*8)/sum(my_subsets(j*8-7:j*8));
    end
    tom_pos(all(isnan(tom_pos),2),:) = [];
    tom_neg(all(isnan(tom_neg),2),:) = [];
    subsets(:,:,i) = [tom_pos tom_neg];
end

subsets = subsets*100;

end