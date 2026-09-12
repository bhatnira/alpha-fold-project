#!/bin/bash
# Submit pipeline as background SLURM job
# Usage: sbatch submit_pipeline.sh

cd /cluster/home/nbhatt04/lean_pipeline

echo "Submitting pipeline job..."
JOB_ID=$(sbatch --parsable run_pipeline.sh)
echo "Job submitted: $JOB_ID"
echo "Monitor with: squeue -j $JOB_ID"
echo "Logs: lean_pipeline_${JOB_ID}.out"
