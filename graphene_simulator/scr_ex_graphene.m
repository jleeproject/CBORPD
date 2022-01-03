chngParamAir = GrapheneChangingParameters(TypeGrapheneSensingSolution.AIR);
chngParamWC0 = GrapheneChangingParameters(TypeGrapheneSensingSolution.WATER_C0);
chngParamWC10 = GrapheneChangingParameters(TypeGrapheneSensingSolution.WATER_C10);
chngParamWC20 = GrapheneChangingParameters(TypeGrapheneSensingSolution.WATER_C20);
fixedParam = GrapheneFixedParameters();


solver = GrapheneModelSolver(fixedParam, chngParamAir);

decVar = GrapheneDecisionVariables();
decVar.Vb = [-80:80];
ids = solver.solve(decVar);

figure(11);clf;plot(ids);

hold on;
ids = solver.setChangingParameters(chngParamWC0).solve(decVar)  ; 
plot(ids);

ids = solver.setChangingParameters(chngParamWC10).solve(decVar)  ; 
plot(ids);

ids = solver.setChangingParameters(chngParamWC20).solve(decVar)  ; 
plot(ids);

sig_vb0 = param_sig_horizontal;    
sig_err = .7e-7;
