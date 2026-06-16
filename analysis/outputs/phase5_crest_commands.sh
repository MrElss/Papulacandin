#!/bin/bash
# CREST conformer ensembles for top Phase-5 candidates (run on QM infra).
# Requires xtb + crest. Papulacandins are highly flexible -> ensembles matter.
# Suggested: GFN2-xTB, implicit water (alpb), to compare exposed hydrophobic SASA.

crest cand01_quinolinecarbonyl.xyz --gfn2 --alpb water --T 8 > cand01.out   # score=0.3958
crest cand02_naphthoyl_6OH.xyz --gfn2 --alpb water --T 8 > cand02.out   # score=0.3847
crest cand03_pyridylphenyl.xyz --gfn2 --alpb water --T 8 > cand03.out   # score=0.3633
crest cand04_biphenyl_4NH2.xyz --gfn2 --alpb water --T 8 > cand04.out   # score=0.3598
crest cand05_biphenyl_4morpholine.xyz --gfn2 --alpb water --T 8 > cand05.out   # score=0.3486
crest cand06_biphenyl_4SO2NH2.xyz --gfn2 --alpb water --T 8 > cand06.out   # score=0.3441
crest cand07_biphenyl_4CONH2.xyz --gfn2 --alpb water --T 8 > cand07.out   # score=0.34
crest cand08_biphenyl_4COOH.xyz --gfn2 --alpb water --T 8 > cand08.out   # score=0.339
crest cand09_biphenyl_4CH2NMe2.xyz --gfn2 --alpb water --T 8 > cand09.out   # score=0.3268
crest cand10_PAPU-0136.xyz --gfn2 --alpb water --T 8 > cand10.out   # score=0.3191
crest cand11_PAPU-0138.xyz --gfn2 --alpb water --T 8 > cand11.out   # score=0.3186
crest cand12_PAPU-0125.xyz --gfn2 --alpb water --T 8 > cand12.out   # score=0.2784
