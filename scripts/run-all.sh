#!/bin/bash
sbatch ./scripts/a/run-n-h.sh
sbatch ./scripts/a/run-o-h.sh
sbatch ./scripts/b/run-n-h.sh
sbatch ./scripts/b/run-o-h.sh
sbatch ./scripts/c/run-n-h.sh
sbatch ./scripts/c/run-o-h.sh
sbatch ./scripts/d/run-n-h.sh
sbatch ./scripts/d/run-o-h.sh
sbatch ./scripts/e/run-n-h.sh
sbatch ./scripts/e/run-o-h.sh
sbatch ./scripts/f/run-n-h.sh
sbatch ./scripts/f/run-o-h.sh

sbatch ./scripts/a/run-n.sh
sbatch ./scripts/a/run-o.sh
sbatch ./scripts/b/run-n.sh
sbatch ./scripts/b/run-o.sh
sbatch ./scripts/c/run-n.sh
sbatch ./scripts/c/run-o.sh
sbatch ./scripts/d/run-n.sh
sbatch ./scripts/d/run-o.sh
sbatch ./scripts/e/run-n.sh
sbatch ./scripts/e/run-o.sh
sbatch ./scripts/f/run-n.sh
sbatch ./scripts/f/run-o.sh
