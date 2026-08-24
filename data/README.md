# Data provenance

`data/raw/Bayesian_polyhyper_muses.hdf5` is downloaded from Zenodo record 13830379 and is intentionally ignored by Git. The required MD5 is `7fd567dccfaea48095ca5df53a8c17d6` and the published size is 60,044,400 bytes.

The compiled production command verifies the checksum before opening HDF5 and selects the largest `log_likelihood` under `posterior_samples/sample*`. For version v1 this is `sample74`; the frozen values are in `configs/pha_map.yaml`.

