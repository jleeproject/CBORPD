function [out] = fn_solve_eq_graphene_LW(Vb, fixedParams, chngParams, constants, L, W, verbose)

    x_init=[0.25 0.3 10^-5];
    
    if(~exist('verbose','var'))
        verbose= false;
    end;
        

    try
        if(verbose)
            yyss=fsolve(@(y)fn_graphene_ext_LW(y, Vb, fixedParams, chngParams, constants, L, W),x_init);
        else
            yyss=fsolve(@(y)fn_graphene_ext_LW(y, Vb, fixedParams, chngParams, constants, L, W),x_init, optimset('Display','off'));
        end
        out = yyss(3);

    catch err
        disp(err)
    end
end