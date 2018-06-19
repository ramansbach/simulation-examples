function [x,fval,exitflag,output] = fmincon_runnested(x0,A,b,Aeq,beq,lb,ub,nonlcon,options,T,dg_vec,T_ref,rho,Lambda_1)

    [x,fval,exitflag,output] = fmincon(@objFun_lnPVec,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    
    function [val,gvals,Hvals] = objFun_lnPVec(x)
        k_max = length(x);
        
       
        ks = 1:k_max;
        val = sum(-(T_ref/T)*(x./ks).*dg_vec +(x./ks).*log((exp(2.5)*ks.^2.5)./(rho*Lambda_1^3*x)));
        val = -val; % since we wish to maximize, and fmincon minimizes
        
        %compute gradient
        gvals = -(T_ref/T)*(1./ks).*dg_vec + (1./ks).*log((exp(2.5)*ks.^2.5)./(rho*Lambda_1^3))-(1./ks).*log(x)-1./ks;
        gvals = -gvals;
        %compute Hessian
        Hvals = zeros(k_max,k_max);
        for k = 1:k_max
           Hvals(k,k) = 1/(k*x(k)); 
        end
        Hvals = -Hvals;
    end
end