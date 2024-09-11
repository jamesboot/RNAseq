sbatch -N 1 -c 1 --partition ncpu --time 72:00:00 -J smrnaseq --wrap "sh /nemo/stp/babs/working/mitterr/projects/patanir/hamish.crerar/hc682_deseq2/smrnaseq/smrnaseq.sh"
