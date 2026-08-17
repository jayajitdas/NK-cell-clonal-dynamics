function plot_fit_with_bootstrap(load_struct)

Ly49H_pos = load_struct.Ly49H_pos;
params = load_struct.k_vals_a0;
model_choice = load_struct.model_choice;
data = collect_subset_data(Ly49H_pos);
params = [params(4) params(3)];

cd27_pos = data(:,1,:)+data(:,2,:)+data(:,5,:)+data(:,6,:);
cd27_neg = data(:,3,:)+data(:,7,:);
mean_cd27_pos = mean(cd27_pos,'all');
mean_cd27_neg = mean(cd27_neg,'all');

switch model_choice
    case {1,2}
        std_pos_data = std(data(:,1:3,:));
        std_neg_data = std(data(:,5:7,:));
        errors = [reshape(std_pos_data,[3 7]); reshape(std_neg_data,[3 7])];
        tom_pos_data = mean(data(:,1:3,:));
        tom_neg_data = mean(data(:,5:7,:));
        data = [reshape(tom_pos_data,[3 7]); reshape(tom_neg_data,[3 7])];
        data = data./sum(data)*100;
    case {5,6,7}
        std_pos_data = [std(data(:,1,:)+data(:,2,:),1) std(data(:,3,:),1)];
        std_neg_data = [std(data(:,5,:)+data(:,6,:),1) std(data(:,7,:),1)];
        errors = [reshape(std_pos_data,[2 7]); reshape(std_neg_data,[2 7])];
        tom_pos_data = [mean(data(:,1,:)+data(:,2,:),1) mean(data(:,3,:),1)];
        tom_neg_data = [mean(data(:,5,:)+data(:,6,:),1) mean(data(:,7,:),1)];
        data = [reshape(tom_pos_data,[2 7]); reshape(tom_neg_data,[2 7])];
        data = data./sum(data)*100;
end

x0 = data(:,1);
t = [3 7 14 23 29 35]-1;

figure()
switch model_choice
    case 1
        sq_params = [1 2 3];
    case 2
        sq_params = [1 2 3];
    case 5
        sq_params = [1 2 3];
    case 6
        sq_params = [1 2 4];
    case 7
        sq_params = [1 2];
end

params = sqrt(params);
fake_errors = ones(length(data(:,1)),length(0:.1:t(end))+1);
fun = @(p,t)call_tamoxifen_model_ODEs_SS(t,p,x0,model_choice,fake_errors,mean_cd27_pos,mean_cd27_neg);
fit = (fun(params,0:.1:t(end)));

for i=1:2
    h(i) = errorbar(1+[0 t],data(i,:),errors(i,:),'LineStyle',"none",'Marker','o','MarkerEdgeColor',"auto",'Linewidth',2,'MarkerSize',10);
    hold on
end
set(gca,'ColorOrderIndex',1)
for i=3:4
    h(i) = errorbar(1+[0 t],data(i,:),errors(i,:),'LineStyle',"None",'Marker','x','MarkerEdgeColor',"auto",'Linewidth',2,'Markersize',10);
    h(i).Bar.LineStyle = 'dashed';
    hold on
end
set(gca,'ColorOrderIndex',1)
plot((0:.1:t(end))+1,fit(1:2,:),'Linewidth',2)

set(gca,'ColorOrderIndex',1)

plot((0:.1:t(end))+1,fit(3:4,:),'--','Linewidth',2)

n_boot = 50000;
for i=1:length(load_struct.predictions_boot(1,1,:))
    for j=1:length(load_struct.predictions_boot(1,:,1))
        uplim(i,j) = prctile(load_struct.predictions_boot(1:n_boot,j,i),97.5);
        lowlim(i,j) = prctile(load_struct.predictions_boot(1:n_boot,j,i),2.5)
    end
end

patch([load_struct.time_fine flip(load_struct.time_fine)]+1,[lowlim(:,1); flip(uplim(:,1))],'b','FaceAlpha',.3)
patch([load_struct.time_fine flip(load_struct.time_fine)]+1,[lowlim(:,2); flip(uplim(:,2))],'r','FaceAlpha',.3)
patch([load_struct.time_fine flip(load_struct.time_fine)]+1,[lowlim(:,3); flip(uplim(:,3))],'b','FaceAlpha',.3)
patch([load_struct.time_fine flip(load_struct.time_fine)]+1,[lowlim(:,4); flip(uplim(:,4))],'r','FaceAlpha',.3)

xlabel('Time (days)','fontweight','bold')
ylabel('%cells','fontweight','bold')
set(gca,'Fontsize',20)
legend('tom+cd27+ data','tom+cd27- data','tom-cd27+ data','tom-cd27- data','tom+cd27+ model','tom+cd27- model','tom-cd27+ model','tom-cd27- model','fontsize',12)

end