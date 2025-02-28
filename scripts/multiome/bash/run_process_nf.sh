#!/bin/bash -l
#SBATCH --output=/scratch/users/%u/%j.out
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH --time=48:00:00

## Runs nfcore rna-seq preprocessing pipeline 
#Arguments:
#   --profile|-p(default:cluster)		Whether to run with a configuration optimised for a slurm HPC cluster (cluster) or local PC (standard)
#   --install|-i(default:true): 		Whether to install nextflow (23.04.1), openjdk (17), and nf-core (if downloadContainers = true)
#   --env|-e;(default:nextflowEnv/):		Name of conda environment to run pipeline, either existing or to be created (if install=true), set to false to run outside an environment
#   "$@"					Add arguments to nextflow run call (e.g. -resume) to end of script call
#Outputs:
#   Pipeline output in processeddata/ directory in cwd 

#Load cuda libraries

#Read the argument values
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -i|--install) install="$2"; shift;;
      -e|--env) env="$2"; shift;;
      -p|--profile) profile="$2"; shift;;
    esac
    shift
done

#Set arguments to default if not specified
install="${install:-true}" #Install packages unless otherwise specified
env="${env:-nextflowEnv/}" #Set env to nextflowEnv/ unless other is specified
profile="${profile:-cluster}" #Set profile to cluster unless other is specified

#Don't install if resume set
if echo "$@" | grep -q "resume"; then
	install="false"
fi

#Install required packages if install=true - Need to create a conda environment to run nfcore on cluster
if [ $install = true ]; then
	conda create --prefix ./${env}
	source ~/.bashrc
	source activate ./${env}
	conda install -y -c bioconda nextflow=23.04.1
	conda install -y -c conda-forge openjdk=17.0.9
	
elif [ ! $env = false ]
	source ~/.bashrc
	source activate ./${env}
fi

nextflow -C scripts/nf/processing/process_nf.config run code/nf/processing/process.nf -profile $profile "$@"

