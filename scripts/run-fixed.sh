#!/bin/bash
sbatch ./scripts/a/run-n.sh
sbatch ./scripts/a/run-n-h.sh
sbatch ./scripts/a/run-o.sh
sbatch ./scripts/a/run-o-h.sh
sbatch ./scripts/b/run-o.sh
sbatch ./scripts/b/run-o-h.sh
sbatch ./scripts/d/run-o.sh
sbatch ./scripts/d/run-o-h.sh
