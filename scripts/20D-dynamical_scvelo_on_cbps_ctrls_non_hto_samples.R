#run scvelo in dynamical mode in cbps hto merged samples
dir.create("outputs/20-RNA_velocity/cbps_ctrl_non_hto")
reticulate::py_run_file("scripts/20D-dynamical_scvelo_on_cbps_ctrl_non_hto_samples.py")

message("Success !")
