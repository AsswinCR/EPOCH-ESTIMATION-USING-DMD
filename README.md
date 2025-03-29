# Epoch Estimation Using Dynamic Mode Decomposition (DMD)

## Project Overview

This repository contains the implementation of a novel approach for Epoch Estimation (also known as Glottal Closure Instant detection) in speech signals using Dynamic Mode Decomposition (DMD). The project was developed for the course "21AIE315: AI in Speech Processing".


## Introduction

During the production of voiced sounds, the lateral vibration of the vocal folds modulates air from the lungs, generating the carrier signal of speech. The significant excitation that occurs in the vocal tract during this process is known as an **Epoch** or **Glottal Closure Instance (GCI)**. Accurate estimation of epoch locations is crucial for determining the correct pitch of speech signals, which serves as an important metric in various speech applications.

Traditional epoch detection algorithms often show degraded performance due to the varying nature of excitation characteristics in speech signals. A key challenge is the interaction between vocal tract response and the excitation signal, which leads to the influence of various vocal tract frequencies mixed with the excitation signal.

## Novelty of Our Approach

Most existing epoch extraction methods assume the speech signal comes from a linear model. However, this assumption is not appropriate for speech signals produced by non-linear interactions between excitation signals and vocal tract responses. Our approach using DMD assumes the speech signal comes from a non-linear speech production system, enabling more accurate epoch location estimation.

Compared to other non-linear methods like Variable Mode Decomposition (VMD), our DMD-based approach has the advantage of using a single decomposition rather than an iterative method, making it computationally more efficient.

## Dataset

We used the CMU Arctic dataset for our experiments, which consists of:
- 1150 utterances from both male and female speakers in US English
- Both speech signals and EGG (Electroglottography) signals
- Recordings in 16-bit, 32kHz sampling frequency

## Methodology

Our approach follows a two-phase framework:

### Phase 1: Finding Instantaneous Pitch Frequencies

1. **Preprocessing**: The speech signal is divided into non-overlapping frames (2000 samples each).
2. **Time-shifting**: The 1D signal vector is converted into a 2D matrix by time-shifting (shift value of 200).
3. **DMD Application**: Dynamic Mode Decomposition is applied to the time-shifted matrix with a rank value of 6.
4. **Mode Frequency Analysis**: 
   - We analyze the speech signal's linear spectrum using Fourier transform to identify sub-signal components.
   - We compare frequencies obtained through DMD with those from the linear spectrum.
   - We select the DMD mode frequency nearest to the average pitch of the speech signal.

### Phase 2: Epoch Estimation

1. **Sine Wave Construction**: Based on the obtained mode frequency, we construct a sine wave that closely resembles the EGG signal.
2. **Epoch Detection**: We identify negative-to-positive zero-crossings in the constructed sine wave as epoch locations.
3. **Validation**: We compare our estimated epochs with the ground truth from the EGG signal to evaluate performance.

## Technical Implementation

### Dynamic Mode Decomposition (DMD)

DMD is a powerful matrix decomposition technique originally developed in fluid dynamics to analyze complex, nonlinear systems without necessarily knowing the underlying governing equations. The procedure involves:

1. **Data Organization**: Arranging time-evolving data snapshots into matrices X1 and X2, where X2 is X1 advanced by one time step.
2. **SVD Decomposition**: Applying Singular Value Decomposition (SVD) to X1.
3. **Low-Rank Approximation**: Computing a low-rank subspace matrix.
4. **Eigendecomposition**: Extracting eigenvalues and eigenvectors to identify dominant modes.
5. **Mode Analysis**: Analyzing the eigenvalues to determine frequency characteristics of each mode.

### Algorithm Flow

1. **Speech Frame Selection**: Extract a frame of 2000 samples from the speech signal.
2. **Pitch Estimation**: Calculate the average fundamental frequency (F0) of the frame.
3. **Spectral Analysis**: Identify dominant frequency peaks from the signal's spectrum.
4. **Window Creation**: Create time-shifted matrices with a specified window length.
5. **DMD Application**: Apply DMD with rank truncation (r=6) to obtain modes.
6. **Mode Selection**: Identify the mode with frequency closest to the average pitch.
7. **Sine Wave Generation**: Construct a sine wave based on the selected mode frequency.
8. **Epoch Identification**: Detect negative-to-positive zero-crossings in the sine wave as epochs.
9. **Concatenation**: Repeat the process for all frames and concatenate results.

## Results

Our experimental results demonstrate that the DMD approach is highly effective in identifying epoch locations accurately when compared to existing methods. The sine wave constructed from the selected mode frequency closely resembles the differentiated EGG signal (DEGG), with negative-to-positive zero-crossings corresponding to actual epoch locations.

Key observations:
- Low-frequency regions in the constructed sine wave correspond to voiced regions where epoch estimation is performed.
- High-frequency regions correspond to unvoiced regions.
- The negative peaks in the DEGG signal closely align with the negative-to-positive zero-crossings in our constructed sine wave.

## Code Structure

The implementation consists of the following components:

1. **Main Script**: Handles the overall processing flow, including signal reading, frame division, and result concatenation.
2. **DMD Functions**: 
   - `DMD_v2`: Creates the time-shifted matrices
   - `dmd_rom`: Performs the core DMD decomposition and returns eigenvalues, eigenvectors, and mode frequencies
3. **Mode Analysis Functions**: 
   - `get_all_ranks`: Extracts rank information across different window lengths
   - `get_freq`: Identifies mode frequencies that match spectral peaks
4. **Signal Generation Functions**:
   - `find_sine`: Generates sine waves based on mode frequencies
   - `simulate_signal`: Orchestrates the complete process for a single frame

## Conclusion

Our work demonstrates that applying the DMD algorithm iteratively on speech signals helps capture the required center frequency of the mode, which is found to be near the fundamental frequency of the glottal excitation signal. This characteristic enables accurate estimation of epoch locations.

Future work could extend this approach to epoch estimation in emotional speech signals, potentially solving problems of low accuracy predictions from existing methods in highly aroused emotions such as anger and happiness.

## References

The project draws upon research in speech signal analysis, linear prediction, phase slope analysis, zero-frequency filtering, and various decomposition methods. Key references include work on DMD by Peter Schmid and various epoch detection techniques like YAGA, SEDREAMS, MMF, and GEFBA.

## Usage

To run the code:

```matlab
% Set the path to your speech data
path = "path_to_your_data/";

% Read the audio file
[y1, Fs] = audioread(path+"your_speech_file.wav");

% Run the epoch detection
all_freqs = simulate_signal(y1(:,1), y1(:,2), start_sample, end_sample, Fs);

% Generate sine wave and plot results
% (See the main script for implementation details)
```

## Dependencies

The code requires MATLAB with the following toolboxes:
- Signal Processing Toolbox
- Audio Toolbox

## License

This project is provided for academic and research purposes. Please cite our work if you use it in your research.
