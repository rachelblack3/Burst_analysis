# Burst_analysis

This code performs various FFT methods upon the raw wavefrom 'burst' data from the Van Allen probes EMFISIS instrument. 

# Theory 

It should be noted early that all technqiues described below to improve our understadning of the raw signals are subject to the uncertainty relation:

$$ \Delta t/\Delta f = 1 $$

## 1. Windowing

Multiplying a signal by a window function in the time domain will result in a convolution of the signal in the frequency domain.

A convolution is an integral operation that gives the inetgral of the product of the two functions, where one is reflected around the y axis and shifted: 

$$ (f*g)(t) = \int _{-\infty }^{\infty }f(\tau )g(t-\tau ) d\tau $$

This is where spectral leakage comes from because the **true signal** is being convolved with the **window**. For example, in the case of a monochromatic signal at X Hz being multipled by a rectangular window with a sample size that doesn't match a multiple of the signal period, the signal will leak across the entire set of Fourier frequencies for the window. 

