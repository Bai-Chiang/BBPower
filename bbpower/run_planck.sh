# #!/bin/bash

sacc_file=/pscratch/sd/k/kwolz/bbdev/SOOPERCOOL/outputs_planck/sacc_files/cl_and_cov_sacc.fits
outdir=output_planck/full_b_diag_cov
mkdir -p $outdir

# Generate fake data
python ./examples/generate_fiducial_spectra.py \
    --globals paramfiles/paramfile_planck.yml \
    --outdir $outdir



# Run pipeline
python -m bbpower BBCompSep --cells_coadded=$sacc_file \
                            --cells_noise="${outdir}/cls_noise.fits" \
                            --cells_fiducial="${outdir}/cls_fid.fits" \
                            --cells_coadded_cov=$sacc_file \
                            --output_dir=$outdir \
                            --config_copy="${outdir}/config.yml" \
                            --config=paramfiles/paramfile_planck.yml

python -m bbpower BBPlotter --cells_coadded_total=$sacc_file \
                            --cells_coadded=$sacc_file \
                            --cells_noise="${outdir}/cls_noise.fits" \
                            --cells_null=$sacc_file \
                            --cells_fiducial="${outdir}/cls_fid.fits" \
                            --param_chains="${outdir}/emcee.npz" \
                            --config=paramfiles/paramfile_planck.yml \
                            --plots="${outdir}/plots"
