ph = discrete_phs(10,0,1)
freqph = frequencies(ph)
set_constraint!(ph, ph.B, ph.D)
ph_noconst  = constraint_elimination(ph)
freqph_no_const = frequencies(ph_noconst)

@test  norm(freqph-freqph_no_const)<1e-10