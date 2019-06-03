# wavelet-denoising
A Matlab implementation for Yu and Guizou method for speech enhancement based on adaptive wavelet denoising on multitaper spectrum.
Final Project of Advanced Methods for Information Representation course, University of Brescia, 2018.

It is well known that the “musical noise” encountered in most frequency domain speech enhancement algorithms is partially due to the large variance estimates of the spectra. To address this issue, we propose in this paper the use of low-variance spectral estimators based on wavelet thresholding the multitaper spectra for speech enhancement. A short-time spectral amplitude estimator is derived which incorporates the wavelet-thresholded multitaper spectra. Listening tests showed that the use of multitaper spectrum estimation combined with wavelet thresholding suppressed the musical noise and yielded better quality than the subspace and MMSE algorithms.

# Usage
Execute SpeechEnhancement.m. You will need to manually change the audio path inside the cose.
To iterate the algorithm:
1) Execute all code once;
2) Repeat sections in the following order:
	A) PSD estimation w/ Multitapers
	B) DWT + Thresholding + Wiener Filter
	C) Reconstruction and results
