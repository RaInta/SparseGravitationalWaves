# SparseGravitationalWaves
Sparse methods (and compressed sensing) applied to gravitational wave signal processing 

### How do I run the code?

Let's see what the main functions do:

```
>> help OMP
   [x, R] = OMP(phi, y, epsilon)
 
  Orthogonal Matching Pursuit recursive reconstruction algorithm
  for compressive sampling.
 
  This solves the inverse problem y = phi*x, by finding x (N X 1 vector),
  where A = M X N measurement matrix, y = M X 1 measurement vector.
  epsilon is a scalar representing the stochastic/residual noise.
  This can be used to define a confidence limit to the Chi-squared distribution.
 
  Very similar to CLEAN routine used in radio interferometry; 
  the stopping criterion is when the residual vector R is reduced below the 
  threshold, epsilon
 

>> help GenSparseProblem
  Initialize random number generator

>> help GenSparseVectors
  GenSparseVectors.m
  [phi, y, s] = GenSparseVectors(M, N, S)
 
  For testing of sparse algorithms.
  This function produces S non-zero coefficients in an N dimensional
  coefficient space from M measurements. 
  Function output is an N-length vector s, with S randomly assigned
  non-zero entries. It also produces a random sensing matrix phi,
  which is created from M measurements of phi*s.
```

Let's generate a measurement vector, y, with 10 measurements, that was the result of a measurement system characterized by a random 10x1024 measurement matrix, phi, from an underlying sparse vector, s, with up to 10 non-zero coefficients:

```
>> [phi, y, s] = GenSparseVectors(10, 1024, 10);
>> plot(1:length(y),y)
>> plot(1:length(s),s)
```

Use OMP to obtain an estimate, x, of the initial vector:

```
>> x = OMP(phi, y);
>> Warning: setting residual error to default [1E-5]
>> length(x)

ans =

   754

>> length(s)

ans =

        1024


```

Note that this version of OMP terminates when it meets the noise conditions imposed upon it, so the length of x here is 754, less than the full vector of coefficients (1024).

Let's look at generating a sparse system, with its own measurement matrix, using `GenSparseVectors`:

```
>> help GenSparseVectors
  GenSparseVectors.m
  [phi, y, s] = GenSparseVectors(M, N, S)
 
  For testing of sparse algorithms.
  This function produces S non-zero coefficients in an N dimensional
  coefficient space from M measurements. 
  Function output is an N-length vector s, with S randomly assigned
  non-zero entries. It also produces a random sensing matrix phi,
  which is created from M measurements of phi*s.

>> [phi, y, s] = GenSparseVectors(10, 1024, 10);
>> x=OMP(phi, y, 1e-6);
>> size(phi)

ans =

         100        1024

>> size(y)

ans =

   100     4

>> size(s)

ans =

        1024           4

>> figure(2)
>> plot(1:length(s),s,'-ko',1:length(x),x,'-rx')
```

### Working with a 4-dimensional wavelet ('Chirplet') bank

Code, similar to that used above, but designed for a 4-dimensional wavelet bank of chirped wave-forms (a 'chirplet') can be run similarly to the above code. The syntax of the functions is the same, except the names are appended with a '4' suffix.

Let's have a look at the chirplets themselves:

```
>> t=0:0.0001:10;
>> y1 =chirpxform(t, f0(1), tau0(1), Q0(1), d0(1));
>> y = zeros(size(t));
    for i = 1:S;
        y = y + chirpxform(t, f0(i), tau0(i), Q0(i), d0(i));
    end
>> plot(t,real(y1))
>> fs=1/(t(2)-t(1))

fs =

       10000

>> specgram(y1,2^4,fs)
>> plot(t,real(y))
>> specgram(y,2^4,fs)
>> xlabel('Time (sec)')
>> ylabel('Amplitude (arbitrary)')
>> title('Chirplet: single injection')
>> gtext({'f_0 = 335 Hz' , '\tau_0 = 3.36 s' , 'Q_0 = 869' , 'd_0 = -2,617 Hz/s'})
```

```
>> help GenSparseVectors4
  GenSparseVectors4.m
  [phi, y, s] = GenSparseVectors4(M, N, S)
 
  For testing of sparse algorithms.
  This function produces S non-zero coefficients in an N dimensional
  coefficient space from M measurements in 4 dimensions.
 
  Function output is an N X 4 matrix s, with S randomly assigned
  non-zero entries. It also produces a random sensing array, phi,
  which is created from M measurements of phi*s.
 
>> plot(1:length(s1),s1,'-ko',1:length(x1),x1,'kx',1:length(s2),s2,'-ro',1:length(x2),x2,'rx',1:length(s3),s3,'-bo',1:length(x3),x3,'bx',1:length(s4),s4,'-go',1:length(x4),x4,'gx')
>> legend('f_0 (inj.)','  (recon.)','\tau','','Q','','d_0','')
>> axis('tight')
>> axis([0 1024 -2.7 1.7])
>> xlabel('Coefficients')
>> ylabel('Amplitude (arbitrary)')
>> title('Four dimensional OMP reconstruction (M_i=100, N_i=1024, S_i=10)')
```
