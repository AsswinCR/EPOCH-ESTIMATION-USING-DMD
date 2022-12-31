
clc;
clear all;
close all;

% Time domain representation of signal
[y1, Fs] = audioread("C:\Users\LENOVO\Desktop\AI2.0\Sem-6\AI in Speech Processing\Research Papers\cmu_us_bdl_arctic\orig\arctic_a0001.wav"); %reading hte wav file
y1=resample(y1,Fs,16000);
sp=y1(:,1);
egg=y1(:,2);
sp = sp(30001:38000);
egg = egg(30001:38000);

%Display
avg_F0 = pitch(sp,Fs);
avg_F0 = avg_F0(1)

X1=[];
Xa=[];
Xb=[];
%window_length=960;
window_length=30;
for i=1:window_length+1
    X1 =[X1;reshape(sp(i:i+length(sp)-window_length-1),[1,length(sp)-window_length])];
end
Xa=X1(1:window_length,:);
Xb=X1(2:window_length+1,:);

number_of_rows = 3000;
X_data = hankel(sp(1:number_of_rows),sp(number_of_rows:end));
Xa=X_data(1:end-1,:);
Xb=X_data(2:end,:);

%[Phi,omega,lambda,b,Xdmd] = DMD(Xa,Xb,8,0.001);

% 1st Decomposition:
Big_X1 = X1;
rank = 4;
selected_mode = 3;
[Big_X1, Phi, ModeFrequencies] = DMD_once(Big_X1, rank, selected_mode, window_length);
plot_signal(rank, avg_F0, sp, egg, Phi, ModeFrequencies);

% 1st Decomposition:
rank = 8;
selected_mode = 7;
[Big_X1, Phi, ModeFrequencies] = DMD_once(Big_X1, rank, selected_mode, window_length);
plot_signal(rank, avg_F0, sp, egg, Phi, ModeFrequencies);

%-------------------------------------------------------------------------

mode_number = 7;
final_mode = real(Phi(mode_number,:));
% Peak estimation for final mode
[pn ln] = findpeaks(final_mode);
ln2 = [];
for i=1:length(ln)
   if(pn(i)>0.04)
       ln2 = [ln2,ln(i)];
   end
end
ln = ln2;
% Peak estimation for DEGG
degg = diff(egg);
[degg_pn degg_ln] = findpeaks(-degg);
degg_ln2 = [];
for i=1:length(degg_ln)
   if(degg_pn(i)>0.025)
       degg_ln2 = [degg_ln2,degg_ln(i)];
   end
end
degg_ln = degg_ln2;

%Accessing the channel 1
figure;
a1=subplot(3,1,1);
plot(sp);axis tight;
xlabel('Time (samples','FontSize',10,'FontWeight','bold');
ylabel('Amplitude','FontSize',10,'FontWeight','bold');
title("Speech signal"+"  ;  "+num2str(avg_F0)+"Hz")

%Accessing the channel 2
a2=subplot(3,1,2);
plot(degg_ln,degg(degg_ln),'r*')
hold on;
plot(diff(egg));axis tight;
xlabel('Time (samples','FontSize',10,'FontWeight','bold');
ylabel('Amplitude','FontSize',10,'FontWeight','bold');
title("DEGG signal")
legend("Original epoch location")

subplot(3,1,3);
epoch = 198:(448-198):2000;
ts = 1:length(epoch);
plot(ln,final_mode(ln),'g*')
hold on;
plot(degg_ln,final_mode(degg_ln),'r*')
hold on;
plot(real(Phi(mode_number,:)));axis tight;
xlabel('Time (samples','FontSize',10,'FontWeight','bold');
ylabel('Amplitude','FontSize',10,'FontWeight','bold');
xlim([0,length(sp)]);
title("Mode - "+num2str(mode_number)+" ;  "+num2str(ModeFrequencies(mode_number))+" Hz");
legend("Max Peaks","Original epoch location");

%-------------------------------------------------------------------------

function [Big_X1, Phi, ModeFrequencies] = DMD_once(Big_X1, rank, selected_mode, window_length)
decomposition_count = 1;
for i = 1:decomposition_count
%rank = 8;
[Eigenvalues, Eigenvectors, ModeAmplitudes, ModeFrequencies, GrowthRates, POD_Mode_Energies]=dmd_rom(Big_X1, rank,1/16000);

Big_X = real(Eigenvectors(selected_mode,:)'); % Select 9th mode for next decomposition
size(Big_X)
Big_X1 = [];
%window_length = round(window_length/2);
for j=1:window_length+1
    Big_X1 =[Big_X1;reshape(Big_X(j:j+length(Big_X)-window_length-1),[1,length(Big_X)-window_length])];
end

disp("Decomposition Number:"+num2str(i))
ModeFrequencies

Phi = Eigenvectors;
end
end

%plot_signal(rank, avg_F0, sp, egg, Phi, ModeFrequencies)
function plot_signal(rank, avg_F0, sp, egg, Phi, ModeFrequencies)
figure;

%Accessing the channel 1
a1=subplot(rank/2+2,1,1);
plot(sp);axis tight;
xlabel('Time (samples','FontSize',10,'FontWeight','bold');
ylabel('Amplitude','FontSize',10,'FontWeight','bold');
title("Speech signal"+"  ;  "+num2str(avg_F0)+"Hz")

%Accessing the channel 2
a2=subplot(rank/2+2,1,2);
plot(diff(egg));axis tight;
xlabel('Time (samples','FontSize',10,'FontWeight','bold');
ylabel('Amplitude','FontSize',10,'FontWeight','bold');
title("DEGG signal")

mode_number = 1;
for i=3:3+(rank/2)-1
    subplot(rank/2+2,1,i);
    plot(real(Phi(mode_number,:)));axis tight;
    epoch = 198:(448-198):2000;
    ts = 1:length(epoch);
    hold on;
    plot(epoch,0.02,'r--*');
    xlabel('Time (samples','FontSize',10,'FontWeight','bold');
    ylabel('Amplitude','FontSize',10,'FontWeight','bold');
    xlim([0,length(sp)]);
    title("Mode - "+num2str(mode_number)+" ;  "+num2str(ModeFrequencies(mode_number))+" Hz");
    mode_number = mode_number+2;
end

end

function plot_signal2(rank, avg_F0, sp, egg, Phi, ModeFrequencies)
figure;

%Accessing the channel 1
a1=subplot(rank+2,1,1);
plot(sp);axis tight;
xlabel('Time (samples','FontSize',10,'FontWeight','bold');
ylabel('Amplitude','FontSize',10,'FontWeight','bold');
title("Speech signal"+"  ;  "+num2str(avg_F0)+"Hz")

%Accessing the channel 2
a2=subplot(rank+2,1,2);
plot(diff(egg));axis tight;
xlabel('Time (samples','FontSize',10,'FontWeight','bold');
ylabel('Amplitude','FontSize',10,'FontWeight','bold');
title("DEGG signal")

mode_number = 1;
for i=3:3+(rank)-1
    subplot(rank+2,1,i);
    plot(real(Phi(mode_number,:)));axis tight;
    xlabel('Time (samples','FontSize',10,'FontWeight','bold');
    ylabel('Amplitude','FontSize',10,'FontWeight','bold');
    xlim([0,length(sp)]);
    title("Mode - "+num2str(mode_number)+" ;  "+num2str(ModeFrequencies(mode_number))+" Hz");
    mode_number = mode_number+1;
end

end


function [Eigenvalues, Eigenvectors, ModeAmplitudes, ModeFrequencies, GrowthRates, POD_Mode_Energies]=dmd_rom(Big_X, r, dt)
    dims=size(Big_X);
    newDims=dims;
    newDims(1)=r;
    %Removes mean. Note: Not removing the mean biases the modes as the
    %data points centroid is shifted. If one wants to capture only the
    %oscillations around the mean, the mean MUST be removed.
    Big_X=Big_X-repmat(mean(Big_X,1),[dims(1) ones(1,length(dims)-1)]);
    
    %Reshapes Big_X
    Big_X=(reshape(Big_X,dims(1),prod(dims(2:end)))).';
    
    
    %Split Big_X into two snapshot sets
    X=Big_X(:,1:end-1);
    Y=Big_X(:,2:end);
        
    %SVD on X
    [U, S, V]=svd(X,0);
    
    %Before reducing rank returns the mode energies for further analysis of
    %the ROM validity
    POD_Mode_Energies=diag(S).^2;
    
    %Reduce rank
    U=U(:,1:r);
    V=V(:,1:r);
    S=S(1:r,1:r);
    
    %Gets A_tilde
    A_tilde=U'*Y*V/S;
    
    %(For debugging), we can compare if A_tilde=A, for r=max(r):
    % A=Y*pinv(X);
    
    %Compute A_tilde eigenvalues and eigenvectors
    [eVecs, Eigenvalues] = eig(A_tilde);
    
    %Gets the DMD eigenvectors back
    Eigenvectors=Y*V*inv(S)*eVecs;
    Eigenvalues=diag(Eigenvalues);
    
    %Gets the mode amplitudes
    ModeAmplitudes=Eigenvectors\X(:,1);
    
    %Gets the frequencies associated with the modes
    fNY=1/(2*dt);
    ModeFrequencies=(angle(Eigenvalues)/pi)*fNY;
    
    %Gets the growth rates    
    GrowthRates=log(abs(Eigenvalues))/dt;
    
    %Reshapes the Eigenvectors back to original Big_X dims
    Eigenvectors=Eigenvectors.';
    Eigenvectors=reshape(Eigenvectors,newDims);
    
end


function centerf = center_freq(sp)
    Fs =32000;
    t = abs(fft(sp));
    t(1) = [];
    tmax = max(t);
    t(ceil(end/2):end) = [];
    abovecutoff = t > tmax / 2;   %3 dB is factor of 2
    lowbin  = find(abovecutoff, 1, 'first');
    highbin = sum(abovecutoff);
    centbin = sqrt(lowbin * highbin);   %geometric mean
    %centbin = mean(lowbin, highbin)
    fft_bin_number = floor(centbin); %FFT bin containing the center frequency
    fft_sp = fft(sp);
    bin_length = round(Fs/length(fft_sp));
    %centerf = fft_bin_number*Fs/length(fft_sp);
    centerf = ((fft_bin_number*Fs/length(fft_sp)) + ((fft_bin_number-1)*Fs/length(fft_sp)))/2;
end