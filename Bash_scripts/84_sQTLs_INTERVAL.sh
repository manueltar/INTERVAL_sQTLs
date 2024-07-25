#!/bin/bash

MASTER_ROUTE=$1
analysis=$2


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0


bashrc_file=$(echo "/home/manuel.tardaguila/.bashrc")

source $bashrc_file
eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""/""$analysis""/")

#rm -rf $output_dir
#mkdir -p $output_dir

Log_files=$(echo "$output_dir""Log_files""/")

rm -rf $Log_files
mkdir -p $Log_files


#### Print_input_file ###################################


type=$(echo "Print_input_file""_""$analysis")
outfile_Print_input_file=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_Print_input_file
echo -n "" > $outfile_Print_input_file
name_Print_input_file=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_Print_input_file=$(echo "$Rscripts_path""290_prepare_input_file.R")


input_file=$(echo "/group/soranzo/manuel.tardaguila/Paper_bits/INTERVAL_sQTLs/Alex_files/sqtl_atu_overlap.tsv")

myjobid_Print_input_file=$(sbatch --job-name=$name_Print_input_file --output=$outfile_Print_input_file --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Print_input_file --input_file $input_file --type $type --out $output_dir")
myjobid_seff_Print_input_file=$(sbatch --dependency=afterany:$myjobid_Print_input_file --open-mode=append --output=$outfile_Print_input_file --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Print_input_file >> $outfile_Print_input_file")

#### split_ceiling ###################################


type=$(echo "split_ceiling""_""$analysis")
outfile_split_ceiling=$(echo "$Log_files""outfile_1.5_""$type"".log")
touch $outfile_split_ceiling
echo -n "" > $outfile_split_ceiling
name_split_ceiling=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_split_ceiling=$(echo "$Rscripts_path""294_split_ceiling_INTERVAL_variants_v2.R")


input_file=$(echo "$output_dir""INTERVAL_sQTL.tsv")
covered_variants=$(echo "$output_dir""covered_variants_DEF.txt")
# --dependency=afterany:$myjobid_Print_input_file

myjobid_split_ceiling=$(sbatch --dependency=afterany:$myjobid_Print_input_file --job-name=$name_split_ceiling --output=$outfile_split_ceiling --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_split_ceiling --input_file $input_file --covered_variants $covered_variants --type $type --out $output_dir")
myjobid_seff_split_ceiling=$(sbatch --dependency=afterany:$myjobid_split_ceiling --open-mode=append --output=$outfile_split_ceiling --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_split_ceiling >> $outfile_split_ceiling")

##### run_LDlinkR  #############################

conda activate renv_latest

type=$(echo "run_LDlinkR_chunk_10")
outfile_run_LDlinkR=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_run_LDlinkR
echo -n "" > $outfile_run_LDlinkR
name_run_LDlinkR=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_run_LDlinkR=$(echo "$Rscripts_path""295_run_LDlinkR_INTERVAL_variants_v2_after_ceiling_split.R")

input_chunk=$(echo "$output_dir""List_chunk_10.txt")
covered_variants=$(echo "/group/soranzo/manuel.tardaguila/Paper_bits/INTERVAL_sQTLs/Run_LDR/covered_variants_DEF.txt")

# # --dependency=afterany:$myjobid_split_ceiling

myjobid_run_LDlinkR=$(sbatch --dependency=afterany:$myjobid_split_ceiling --job-name=$name_run_LDlinkR --output=$outfile_run_LDlinkR --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_run_LDlinkR --input_chunk $input_chunk --covered_variants $covered_variants --type $type --out $output_dir")
myjobid_seff_run_LDlinkR=$(sbatch --dependency=afterany:$myjobid_run_LDlinkR --open-mode=append --output=$outfile_run_LDlinkR --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_run_LDlinkR >> $outfile_run_LDlinkR")

conda deactivate

#### Collect_LDlinkR_results #############################


type=$(echo "Collect_LDlinkR_results")
outfile_Collect_LDlinkR_results=$(echo "$Log_files""outfile_3_""$type"".log")
touch $outfile_Collect_LDlinkR_results
echo -n "" > $outfile_Collect_LDlinkR_results
name_Collect_LDlinkR_results=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_Collect_LDlinkR_results=$(echo "$Rscripts_path""296_Gather_variants_from_LDlinkR_INTERVAL.R")

input_file=$(echo "$output_dir""INTERVAL_sQTL.tsv")
LD_thresholds=$(echo '0.5,0.7,0.9')

# --dependency=afterany:$myjobid_run_LDlinkR

myjobid_Collect_LDlinkR_results=$(sbatch --dependency=afterany:$myjobid_run_LDlinkR --job-name=$name_Collect_LDlinkR_results --output=$outfile_Collect_LDlinkR_results --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Collect_LDlinkR_results --input_file $input_file --LD_thresholds $LD_thresholds  --type $type --out $output_dir")
myjobid_seff_Collect_LDlinkR_results=$(sbatch --dependency=afterany:$myjobid_Collect_LDlinkR_results --open-mode=append --output=$outfile_Collect_LDlinkR_results --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Collect_LDlinkR_results >> $outfile_Collect_LDlinkR_results")
