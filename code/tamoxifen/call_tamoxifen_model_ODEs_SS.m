function output = call_tamoxifen_model_ODEs_SS(time,p,x0,model_choice,errors,I,M)

kM = p(1)^2;
if model_choice == 7
    kI = p(2)^2;
elseif model_choice == 8
    kI = 0;
end
r = kM*M/I;
lambda = kM*M-kI*I;
fun = @(t,x)tamoxifen_model_ODEs(t,x,[lambda,r,kM,kI],model_choice);
sol = ode45(fun,[0 time(end)],x0);

output = (deval(sol,time))./errors(:,2:end);

end