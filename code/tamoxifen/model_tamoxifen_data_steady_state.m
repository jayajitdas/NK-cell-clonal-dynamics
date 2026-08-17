function model_tamoxifen_data_steady_state()

% Created 8/3/23
% Constraining parameters to equalities

n_bootstrap = 1000;
model_choices = [1 2 5 6 1 2 5 6];
model_choices = [5 6 5 6];
model_choices = [7 8];% 7 8];

for z=1:length(model_choices)
    clear k_vals chi2 resid exitflag
    model_choice = model_choices(z);
    %if z<=length(model_choices)/2
    %    Ly49H_pos = true;
    %else
        Ly49H_pos = false;
    %end
    
    raw_data = collect_subset_data(Ly49H_pos);


data = raw_data;

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
        disp('Not yet fixed')
    case {5,6,7,8}
        std_pos_data = [std(data(:,1,:)+data(:,2,:),1) std(data(:,3,:),1)];
        std_neg_data = [std(data(:,5,:)+data(:,6,:),1) std(data(:,7,:),1)];
        errors = [reshape(std_pos_data,[2 7]); reshape(std_neg_data,[2 7])];
        tom_pos_data = [mean(data(:,1,:)+data(:,2,:),1) mean(data(:,3,:),1)];
        tom_neg_data = [mean(data(:,5,:)+data(:,6,:),1) mean(data(:,7,:),1)];
        data = [reshape(tom_pos_data,[2 7]); reshape(tom_neg_data,[2 7])];
        data = data./sum(data)*100;
end

x0 = data(:,1);
y_fit = (data(:,2:end))./errors(:,2:end);
t = [3 7 14 23 29 35]-1;

fun = @(p,t)call_tamoxifen_model_ODEs_SS(t,p,x0,model_choice,errors,mean_cd27_pos,mean_cd27_neg);
switch model_choice
    case 1
        init = 0.1*ones(4,1);
        sq_params = [1 2 3];
    case 2
        init = 0.01*ones(6,1);
        sq_params = [1 2 3];
    case 3
        init = 0.1;
    case 4
        init = 0.1*ones(3,1);
    case 5
        init = 0.1*ones(3,1);
        sq_params = [1 2 3];
    case 6
        init = 0.1*ones(4,1);
        sq_params = [1 2 4];
    case 7
        init = [0.1 0.1];
        sq_params = [];
    case 8
        init = 0.1;
        sq_params = [];
end

solve_options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','FunctionTolerance',1e-10,'MaxFunctionEvaluations',1e8,'MaxIterations',1e6,'Display','none');
[fit_p_a0,chi2_a0,resid,exitflag_a0,~,~,J] = lsqcurvefit(fun,init,t,y_fit,[],[],solve_options);
kM = fit_p_a0(1)^2;
if model_choice == 7
    kI = fit_p_a0(2)^2;
elseif model_choice == 8
    kI = 0;
end
r = kM*mean_cd27_neg/mean_cd27_pos;
lambda = kM*mean_cd27_neg-kI*mean_cd27_pos;
k_vals_a0 = [lambda,r,kI,kM];

for i=1:n_bootstrap
    data = zeros(size(raw_data));
    for j=1:size(raw_data,3)
        ind = randi(size(raw_data,1),[size(raw_data,1) 1]);
        data(:,:,j) = raw_data(ind,:,j);
    end
switch model_choice
    case {1,2}
        std_pos_data = std(data(:,1:3,:));
        std_neg_data = std(data(:,5:7,:));
        errors = [reshape(std_pos_data,[3 7]); reshape(std_neg_data,[3 7])];
        tom_pos_data = mean(data(:,1:3,:));
        tom_neg_data = mean(data(:,5:7,:));
        data = [reshape(tom_pos_data,[3 7]); reshape(tom_neg_data,[3 7])];
        data = data./sum(data)*100;
    case {5,6,7,8}
        std_pos_data = [std(data(:,1,:)+data(:,2,:),1) std(data(:,3,:),1)];
        std_neg_data = [std(data(:,5,:)+data(:,6,:),1) std(data(:,7,:),1)];
        errors = [reshape(std_pos_data,[2 7]); reshape(std_neg_data,[2 7])];
        tom_pos_data = [mean(data(:,1,:)+data(:,2,:),1) mean(data(:,3,:),1)];
        tom_neg_data = [mean(data(:,5,:)+data(:,6,:),1) mean(data(:,7,:),1)];
        data = [reshape(tom_pos_data,[2 7]); reshape(tom_neg_data,[2 7])];
        data = data./sum(data)*100;
end

x0 = data(:,1);
y_fit = (data(:,2:end))./errors(:,2:end);
t = [3 7 14 23 29 35]-1;

fun = @(p,t)call_tamoxifen_model_ODEs_SS(t,p,x0,model_choice,errors,mean_cd27_pos,mean_cd27_neg);

%solve_options = optimoptions('lsqcurvefit','Display','none','MaxFunctionEvaluations',10000);%,'FunctionTolerance',1e-10,'MaxFunctionEvaluations',1e8,'MaxIterations',1e6,'Display','none');
[fit_p,chi2(i),resid,exitflag(i),~,~,J] = lsqcurvefit(fun,init,t,y_fit,[],[],solve_options);

kM = fit_p(1)^2;
if model_choice == 7
    kI = fit_p(2)^2;
elseif model_choice == 8
    kI = 0;
end
r = kM*mean_cd27_neg/mean_cd27_pos;
lambda = kM*mean_cd27_neg-kI*mean_cd27_pos;
k_vals(i,:) = [lambda,r,kI,kM];

    fit_p(sq_params) = sqrt(fit_p(sq_params));

    % Use fake errors to call ODE
    fake_errors = ones(length(x0), n_time);
    fun_pred = @(p,t)call_tamoxifen_model_ODEs(t, p, x0, model_choice, fake_errors);

    predictions_boot(i,:,:) = fun_pred(fit_p, time_fine);
end

k_vals(:,sq_params) = k_vals(:,sq_params).^2;
k_vals_a0(sq_params) = k_vals_a0(sq_params).^2;
boo = any(imag(k_vals)==0,2);
ci{z} = [prctile(real(k_vals(boo,:)),2.5,1)' k_vals_a0' prctile(real(k_vals(boo,:)),97.5,1)'];

end

end