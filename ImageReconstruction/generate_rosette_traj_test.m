kmax = 500;
w1 = 1611; % in rads/second
w2 = w1;
sampling_interval = 10e-6; % in seconds
samples_per_petal = 432;
phi = -pi/4;
beta = 0;
t = (0:samples_per_petal-1) * sampling_interval/2;
for phi=-pi/2:0.1:pi/2
    for beta=0:0.1:pi
        kxy = kmax*cos(phi)*sin(w1*t)*exp(1i*w2+beta);
        kx(ii,jj,kk) = real(kxy);
        ky = imag(kxy);
        kz = kmax*sin(phi)*sin(w1*t);
    end
end
