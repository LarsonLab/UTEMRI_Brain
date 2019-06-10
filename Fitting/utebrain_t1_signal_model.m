function Smodel = utebrain_t1_signal_model(x, num_components, TE, flip, TR)
num_variables = 5;
% x contains M0, T2, df, phi, T1
Smodel = 0;
for n = 0:num_components-1
    Smodel = Smodel + ...
        x(num_variables*n+1) .* ... % M0
        sin(flip) .* (1-exp(-TR/x(num_variables*n+5))) ./ (1-cos(flip)*exp(-TR/x(num_variables*n+5))) .*  ...  %T1 weighting
        exp( -TE(:)/x(num_variables*n+2)-i*2*pi*x(num_variables*n+3)*TE(:) + i*x(num_variables*n+4));  % T2, df weighting
end

end

