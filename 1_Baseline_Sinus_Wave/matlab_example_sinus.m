num = fir1(31,0.5);
fir = dsp.FIRFilter('Numerator',num);  
iir = dsp.IIRFilter('Numerator',sqrt(0.75),...
        'Denominator',[1 -0.5]);
x = iir(sign(randn(2000,25))); 
n = 0.1*randn(size(x));           
d = fir(x) + n; 

l = 20;  % order
mu = 0.0510;  %tuneing parameter
m  = 5;  % the decimation factor for analysis and simulation

figure;
lms = dsp.LMSFilter('Length',l,'StepSize',mu);
[mmse,emse,meanW,mse,traceK] = msepred(lms,x,d,m);
[simmse,meanWsim,Wsim,traceKsim] = msesim(lms,x,d,m);

nn = m:m:size(x,1);
semilogy(nn,simmse,[0 size(x,1)],[(emse+mmse)...
    (emse+mmse)],nn,mse,[0 size(x,1)],[mmse mmse])
title('Mean Squared Error Performance')
axis([0 size(x,1) 0.001 10])
legend('MSE (Sim.)','Final MSE','MSE (Pred.)','Min. MSE')
xlabel('Time Index')
ylabel('Squared Error Value')