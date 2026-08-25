# Literature audit

Audit date: 2026-08-25. Searches covered arXiv, INSPIRE, and publisher pages through August 2026. The table records direct scope overlap, not an exhaustive citation graph. No priority claim is made.

| Reference | Model | Finite density | Critical point | k=0 | k!=0 | Full spectrum | Stable/metastable/spinodal | Numerical method | Direct overlap with this work |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| Hippert et al., PRD 110 (2024), arXiv:2309.00579 | Bayesian PHA/PA EMD | yes | yes | no | no | no | thermodynamics | Bayesian EMD backgrounds | exact model and posterior |
| Khan et al., arXiv:2603.20482 | Bayesian PHA/PA EMD | yes | yes | no | no | no | stable/first-order thermodynamics | native C++ backgrounds/transport | exact backgrounds, UV extraction, external transport |
| Rougemont, Critelli, Noronha, PRD 98 (2018), arXiv:1804.00189 | earlier QCD-calibrated EMD | yes | yes | yes | no | lowest homogeneous modes | selected stable paths | shooting | closest homogeneous benchmark; different EMD realization |
| Alho et al., PRD 101 (2020), arXiv:2002.09544 | improved holographic QCD | no | deconfinement transition | scalar | no | selected scalar modes | stable black holes | shooting/spectral | homogeneous temperature trajectories in a different model |
| DeWolfe, Gubser, Rosen, PRD 84 (2011), arXiv:1108.2029 | holographic critical-point EMD | yes | yes | transport limit | hydrodynamic limit | no | phase diagram | Kubo/linear response | diffusion vanishing and dynamic critical behavior |
| Janik, Jankowski, Soltanpanahi, PRL 117 (2016), arXiv:1512.06871 | Einstein--scalar plasma | no conserved density | first-order transition | yes | yes | low modes | spinodal | QNM boundary-value problem | sound instability and preferred scale |
| Janik, Jankowski, Soltanpanahi, JHEP 06 (2016) 047, arXiv:1603.05950 | Einstein--scalar plasma | no conserved density | phase transitions | yes | yes | low modes | unstable branches | QNM boundary-value problem | hydro/nonhydro dynamics near transitions |
| Kovtun and Starinets, PRD 72 (2005), hep-th/0506184 | AdS5 black brane | no | no | yes | yes | channel spectra | stable | gauge-invariant ODEs | standard source-free QNM prescription |
| Kaminski et al., JHEP 02 (2010) 021, arXiv:0911.3610 | coupled holographic fields/D3-D7 | yes | no | yes | yes | coupled poles | stable | coupled Green functions | operator mixing and coupled-field diagnostics |
| Jansen, EPJ Plus 132 (2017), arXiv:1709.09178 | general black-hole perturbations | model dependent | model dependent | yes | yes | yes | model dependent | pseudospectral generalized eigenproblem | reference numerical architecture |
| Abbasi and Tahery, JHEP 10 (2020) 076, arXiv:2007.10024 | charged AdS5 Reissner--Nordstrom plasma | yes | no | yes | yes | spin-0/1/2 spectra | stable | coupled spectral QNMs/complex momentum | charged pole collisions and hydrodynamic convergence |
| Cruz Rojas, Demircik, Järvinen, PRD 111 (2025) 046017, arXiv:2405.02399 | V-QCD | yes | quantum critical point | yes | yes | coupled QNMs | modulated unstable region | spectral QNMs | broad finite-density spectrum in a different model |
| Villarreal et al., PRD 112 (2025) 066017, arXiv:2507.04165 | Einstein--dilaton holographic QCD | no baryon charge | no CEP | yes | selected finite momentum | tensor/vector/scalar | stable/small black holes | pseudospectral | recent QCD-like tensor/vector comparison |
| Betzios et al., PRD 97 (2018) 081901 | Einstein--dilaton critical plasma | no charge | continuous transition | yes | partly | analytic low spectrum | critical branch | analytic + numerical | accumulation/slow nonhydrodynamic modes in another universality setting |
| Guo, Kuang, Qian, JHEP 06 (2025) 142, arXiv:2410.05065 | EMD critical black holes | yes | yes | perturbative stability | limited | no | underlying branches | thermodynamic/geometric analysis | recent EMD phase-structure overlap, not PHA spectroscopy |

## Audit conclusion

The literature establishes homogeneous QNMs in an earlier QCD-calibrated EMD realization, improved-holographic-QCD temperature trajectories, and broader finite-momentum spectra in other models. The Route-B manuscript distinguishes itself by using the public PHA MAP realization, reconstructing its phase structure, independently validating all reported homogeneous poles, and following all three sectors on controlled thermodynamic trajectories.  As a quantitative external benchmark, the present pipeline also reproduces the nine approximate relaxation times quoted in Rougemont--Critelli--Noronha (2018) and compares both models at two common neutral states and their model-specific endpoints. This is a scope comparison, not a priority claim; the search should be repeated immediately before submission.

## Primary links

- https://arxiv.org/abs/2309.00579
- https://zenodo.org/records/13830379
- https://arxiv.org/abs/2603.20482
- https://arxiv.org/abs/1804.00189
- https://arxiv.org/abs/1108.2029
- https://arxiv.org/abs/1512.06871
- https://arxiv.org/abs/1603.05950
- https://arxiv.org/abs/hep-th/0506184
- https://arxiv.org/abs/0911.3610
- https://arxiv.org/abs/1709.09178
- https://arxiv.org/abs/2007.10024
- https://arxiv.org/abs/2002.09544
- https://arxiv.org/abs/2405.02399
- https://arxiv.org/abs/2507.04165
- https://arxiv.org/abs/2410.05065
