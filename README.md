
## EPOCH ESTIMATION USING DMD

📝 Project Overview

This project aims to estimate epoch or Glottal Closure Instant (GCI) locations from speech signals using Dynamic Mode Decomposition (DMD). Unlike conventional linear models, DMD models the speech production system as a non-linear system to capture instantaneous pitch frequencies from the modes of DMD, providing higher accuracy in epoch detection.



🎯 Objective
	•	Estimate epoch locations with higher accuracy using DMD.
	•	Identify the mode frequency closest to the fundamental frequency for each glottal cycle.
	•	Construct a sine wave from the estimated mode and extract negative-to-positive zero-crossings as epoch locations.
	•	Compare DMD-based epoch estimation with traditional methods such as VMD and ZFF.

⸻

Key Features:
	•	✅ Non-linear Model for Speech Production: DMD captures the non-linear interaction of the excitation signal with vocal tract responses.
	•	✅ Accurate Epoch Estimation: Detects epoch locations with greater accuracy by choosing the dominant mode.
	•	✅ Efficient Pipeline: Single-step DMD decomposition is computationally lighter than VMD’s iterative process.

⸻

📂 File Structure

├── DMD_code_5.m                # MATLAB code for DMD-based epoch estimation
├── README.md                   # Project documentation
├── Report.docx                  # Detailed project report with results and analysis
└── SPEECH_DMD_FINAL.pptx       # Final presentation with methodology and results



⸻

📖 Usage Instructions

1️⃣ Clone the Repository

git clone https://github.com/AsswinCR/EPOCH-ESTIMATION-USING-DMD.git
cd EPOCH-ESTIMATION-USING-DMD



⸻

2️⃣ Run the MATLAB Code
	•	Open DMD_code_5.m in MATLAB.
	•	Execute the script to process the speech signal and estimate epoch locations.
	•	Modify the input parameters (if needed) to experiment with different datasets.

% Run in MATLAB
DMD_code_5.m



⸻

3️⃣ Dataset Description
	•	Dataset: CMU Arctic Dataset
	•	Details:
	•	1150 utterances of male and female speakers.
	•	US English language.
	•	Speech signals and EGG signals.
	•	Recorded at 32KHz sampling frequency with 16-bit resolution.

⸻

⚙️ Installation and Requirements

📚 Required Software
	•	MATLAB (R2021 or later recommended)
	•	Toolboxes:
	•	Signal Processing Toolbox
	•	Optimization Toolbox

⸻

Step-by-Step Implementation

Step 1: Load and Preprocess the Data
	•	Load the speech signal and corresponding EGG signal.
	•	Resample both signals to 16KHz.

[y1, Fs] = audioread(path+"arctic_a0002.wav");
sp1 = y1(:,1); % Speech signal
egg1 = y1(:,2); % EGG signal
sp = resample(sp1, Fs, 16000);
egg = resample(egg1, Fs, 16000);



⸻

Step 2: Apply DMD to Extract Modes
	•	Decompose the speech signal into different modes using DMD.
	•	Identify the mode frequency closest to the fundamental frequency.

% Define parameters
start_limit = 20000;
end_limit = 52000;

% Iterate through speech frames
for k = start_limit:2000:end_limit-2000
    range1 = k;
    range2 = k + 2000;
    final_mode_freq = simulate_signal(sp, egg, range1, range2, Fs);
    all_freqs = [all_freqs final_mode_freq];
end



⸻

Step 3: Generate Sine Wave from Mode Frequency
	•	Construct a sine wave based on the estimated mode frequency.
	•	Combine sine waves from all speech frames to form the simulated excitation signal.

% Generate sine wave using mode frequency
fs = 16000;
dt = 1/fs;
StopTime = ((62.5 * 2 * length(all_freqs)) / 1000);
t = (0:dt:StopTime)';
data_all = [];

for i = 1:length(all_freqs)
    final_mode_freq = all_freqs(i);
    data = find_sine(final_mode_freq);
    data_all = [data_all; data];
end

data_all = data_all(1:length(data_all)-length(all_freqs)+1);



⸻
Step 4: Epoch Estimation
	•	Identify epoch locations from the zero-crossings of the sine wave.
	•	Compare epoch locations with the DEGG signal to validate accuracy.

% Plot DEGG and sine wave to visualize epoch locations
figure;
egg_part = egg(start_limit:end_limit);
plot(diff(egg_part));
hold on;
plot(t * 100000 / 6.25, data_all);



⸻

Step 5: Mode Frequency Estimation using DMD
	•	Apply DMD to speech frames to estimate dominant mode frequencies.
	•	Identify modes that align closely with the fundamental frequency.

function [final_mode_freq] = simulate_signal(sp, egg, range1, range2, Fs)
    sp = sp(range1:range2-1);
    egg = egg(range1:range2-1);

    % Apply DMD and estimate mode frequency
    [X1, X2] = DMD_v2(sp, 200);
    [mode_freq] = get_freq(X1, X2, Fs);
    final_mode_freq = mode_freq;
end



⸻

📊 Results and Discussion

✅ Experimental Setup
	•	Dataset: CMU Arctic Dataset with 1150 utterances of male and female speakers.
	•	Sampling Frequency: 32KHz (downsampled to 16KHz).
	•	Validation: Comparison with ground truth EGG signals.

⸻

📈 Observations
	•	DMD accurately identifies mode frequencies closer to the fundamental frequency.
	•	Estimated epoch locations closely resemble the DEGG signal’s negative peaks.
	•	Computational complexity reduced compared to iterative VMD methods.

⸻

🎥 Visualizations
	•	Epoch locations plotted alongside DEGG signals to validate accuracy.
	•	Multiple speech frames concatenated to visualize complete speech signals.

⸻

🤝 Contribution Guidelines

Contributions, improvements, and suggestions are welcome!
To contribute:
	1.	Fork the repository.
	2.	Create a new branch (feature/your-feature).
	3.	Commit your changes.
	4.	Submit a pull request.

⸻

🧩 License
This project is licensed under the MIT License.

⸻

🗂️ References
	•	Karam, Zahi N., et al. “Ecologically valid long-term mood monitoring using speech.” ICASSP 2014.
	•	Naylor, P.A., et al. “Estimation of glottal closure instants in voiced speech using the DYPSA algorithm.” IEEE 2007.
	•	Dragomiretskiy, K., Zosso, D. “Variational mode decomposition.” IEEE Transactions 2014.


