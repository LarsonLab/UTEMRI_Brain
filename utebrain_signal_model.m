    function Smodel = utebrain_signal_model(x, num_components, TE)
        Smodel = 0;
        for n = 0:num_components-1
            Smodel = Smodel + x(4*n+1) * exp( -TE(:)/x(4*n+2)-i*2*pi*x(4*n+3)*TE(:) + i*x(4*n+4));
        end            
        
    end

