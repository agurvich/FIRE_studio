;-------------------------------------------------------------------------------------
; Pei (1992) MW,LMC,SMC dust laws used to calculate differential extinction
;   for dust at most frequencies
; (INPUT LAMBDA IN MICRONS)
;     NOTE: Pei's model is based on observations from lambda > 1000 Angstroms
;       (f < 3.0e15 Hz), but the grain models extrapolate these fits at least to
;		lambda = 0.001 micrometer (= 10 Angstroms, or f = 3.0e17 Hz), where 
;		the extinction drops off rapidly. So it's probably safe to use for all freq.
;-------------------------------------------------------------------------------------
pro pei_dustparam, lambda, xsi, MW=MW, LMC=LMC, SMC=SMC

if keyword_set(MW) then begin
  a = [165., 14., 0.045, 0.002, 0.002, 0.012]
  l = [0.047, 0.08, 0.22, 9.7, 18., 25.]
  b = [90., 4.00, -1.95, -1.95, -1.80, 0.00]
  n = [2.0, 6.5, 2.0, 2.0, 2.0, 2.0]
  R_V = 3.08
endif
if keyword_set(LMC) then begin
  a = [175., 19., 0.023, 0.005, 0.006, 0.020]
  l = [0.046, 0.08, 0.22, 9.7, 18., 25.]
  b = [90., 5.50, -1.95, -1.95, -1.80, 0.00]
  n = [2.0, 4.5, 2.0, 2.0, 2.0, 2.0]
  R_V = 3.16
endif
if keyword_set(SMC) then begin
  a = [185., 27., 0.005, 0.010, 0.012, 0.030]
  l = [0.042, 0.08, 0.22, 9.7, 18., 25.]
  b = [90., 5.50, -1.95, -1.95, -1.80, 0.00]
  n = [2.0, 4.0, 2.0, 2.0, 2.0, 2.0]
  R_V = 2.93
endif

xsi = 0*lambda
for i = 0, 5 do begin
  xsi = xsi + a[i] / ( (lambda/l[i])^(n[i]) + (l[i]/lambda)^(n[i]) + b[i] )
endfor
R_lam = (1.0 + R_V) * xsi
end
