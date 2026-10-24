if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi
cd /home/ravanox/codex/bioinformatics
source /home/ravanox/codex/bioinformatics/tools/miniforge3/etc/profile.d/conda.sh
export PATH="/home/ravanox/codex/bioinformatics/tools/samtools-env/bin:/home/ravanox/codex/bioinformatics/tools/gatk-4.7.0.0:$PATH"
case "$BIO_TOOL" in
    bwa) conda activate /home/ravanox/codex/bioinformatics/tools/miniforge3/envs/bwa; bwa ;;
    fastqc) conda activate /home/ravanox/codex/bioinformatics/tools/miniforge3/envs/fastqc; fastqc --help ;;
    gatk) gatk --version ;;
    samtools) samtools --version ;;
    miniforge) conda activate base; conda --version ;;
esac
unset BIO_TOOL
