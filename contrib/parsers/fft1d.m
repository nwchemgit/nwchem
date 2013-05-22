##
## fft1d.m
##
## Simple GNU Octave script for doing a 1D forward Fourier Transform.
## Tested for Octave 3.2.4 on GNU/Linux, no guarantees will work for
## matlab...
##
## Example usage (note "-q" for silent so you can redirect output):
##
## octave -q fft1d.m ft.dat > fw.dat
##
## where the input (ft.dat) has data: t, f(t) and the output (fw.dat)
## will have w, Re[f(w)], Im[f(w)], Abs[f(w)].
##
## You will have to manually change the flag, damping, padding in this
## file.
##
## XXXTODO command line switches for opts
##
## Kenneth Lopata
##
## Last modified: May 22, 2013
##


##
## Switches and options for preprocessing
##
preprocess_zero = true;       #zero time signal at t=0 (req'd doing padding or damping)
preprocess_pad  = false;
preprocess_damp = false;

preprocess_expconst = 50.0;   #if damping: damp by exp(-t/tau) before FFT; same time units as input
preprocess_npad = 50000;      #if padding: add this many points to time signal before FFT


##
## Read in raw data t, f(t) from file (1st command line arg)
##
arglist = argv ();
fname = arglist{1};
data = load (fname);
t = data (:,1);
f = data (:,2);

##
## Optional preprocessing of time signal
##
## (zero at t=0)
if preprocess_zero
  f0 = f(1);
  f = f - f0;
endif

## (exponential damping)
if preprocess_damp
  damp = exp ( -(t-t(1)) / preprocess_expconst);
  f = f .* damp;
endif

## (zero padding)
if preprocess_pad
  zeros = linspace (0.0, 0.0, preprocess_npad);
  f = vertcat (f, zeros(:));
endif


##
## Do FFT, compute frequencies, and print to stdout: w, Re(fw), Im(fw),
## Abs(fw). Note we only print the positive frequencies (first half of
## the data)--this is fine since time signal is pure real so the
## negative frequencies are the Hermitian conjugate of the positive
## frequencies.
##
fw = fft (f);

n = size (f,1);     #this includes padding
dt = t(2) - t(1);   #assumes constant time step; XXX no safety checks
period = (n-1)*dt - t(1)
dw = 2.0 * pi / period;

m = n / 2;
wmin = 0.0;
wmax = m*dw;

fw_pos = fw(1:m);              #FFT values of positive frequencies (first half of output array)
fw_re = real(fw_pos);
fw_im = imag(fw_pos);
fw_abs = abs(fw_pos);

w = linspace (wmin, wmax, m);  #positive frequency list

output = horzcat (w(:), fw_re, fw_im, fw_abs);
format short e;
disp (output);
