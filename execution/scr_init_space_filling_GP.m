initializer = BoFactory.getInitializer(TypeInitializer.Sobol, simulator, xDomain, n_init_sample, setting.optProb.nCon, setting.samplesize_for_variance_init);
fprintf('Drawing initial (%d x %d) samples.....', initializer.samplesizeForVariance, initializer.nInitSample);
[yy_arr, xx_arr, samplesize_arr, yy_cons_arr, y_orig_cell] = initializer.getSamples();
fprintf('Done.\n')


