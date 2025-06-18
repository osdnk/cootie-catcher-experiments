#!/bin/bash
sbatch ./scripts/c/run-n-h.sh
sbatch ./scripts/c/run-o-h.sh
sbatch ./scripts/d/run-n-h.sh
sbatch ./scripts/d/run-o-h.sh
sbatch ./scripts/c/run-n.sh
sbatch ./scripts/c/run-o.sh
sbatch ./scripts/d/run-n.sh
sbatch ./scripts/d/run-o.sh
