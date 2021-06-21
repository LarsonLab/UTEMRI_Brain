function Smodel = UTE_T1T2_signal_model(x, TE, flip, general_opts)

num_components = general_opts.num_components;
TR = general_opts.TR;

num_variables = 5;
% x contains M0, T2, df, phi, T1
Smodel = 0;
for n = 0:num_components-1
    M0 = x(num_variables*n+1);
    T2 = x(num_variables*n+2);
    df = x(num_variables*n+3);
    if n == 1 && general_opts.fixDf
        df = x(num_variables*(n-1)+3)+general_opts.methylene_freq_est;
    end
    phi = x(num_variables*n+4);
    T1 = x(num_variables*n+5);
    Smodel = Smodel + ...
        M0 .* ... % M0
        sin(flip) .* (1-exp(-TR/T1)) ./ (1-cos(flip)*exp(-general_opts.TR/T1)) .*  ...  %T1 weighting
        exp(-TE(:)/T2-1i*2*pi*df*TE(:) + 1i*phi);  % T2, df weighting
end

end

