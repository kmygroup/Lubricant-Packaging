function [x,fval,exitflag,output,lambda,grad,hessian] = opt_fmin_solver_v2(x0,lb,ub)
x0=[0.005*ones(291,1);]; %initial guess for decision variables

lb=[0*ones(291,1);];
ub=[0.00303*ones(291,1);];

% Start with the default options
options = optimset;
% Modify options setting
options = optimset(options,'Display' ,'iter');
options = optimset(options,'PlotFcns' ,{  @optimplotx @optimplotfval });
options = optimset(options,'Diagnostics' ,'on');
options = optimset(options,'LargeScale' ,'on');
options = optimset(options,'TolX',1e-35,'TolFun',1e-35,'TolCon',1e-35);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@mult_obj_optflush_v2,x0,[],[],[],[],lb,ub,[],options);
save output_mult_obj_optflush
end


