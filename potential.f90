function v(x)
!analitical function of the potential (Morse Potential)
!the units of the potential must be in J.

implicit none

!De: well depth
!a: control of the width of the potential 
!re: equilibrium distance
!x: positions of the oscillator

real(8) :: De, a, re, v,x

re=1.275d-10 !m
De=7.10647d-19 !J
a=1.81181d10

v=De*(1-exp(-a*(x-re)))**2
endfunction
