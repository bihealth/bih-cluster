# Slurm Rosetta Stone

!!! faq "Rosetta Stone?"
    The [Rosetta Stone](https://en.wikipedia.org/wiki/Rosetta_Stone) is a stone slab that carries the same text in Egyptian hieroglyphs and ancient Greek.
    This was key for decyphering Egyptian hieroglyphs in the 18th century.
    Nowadays, the term is often used to label translation tables such as the one below.

The table below shows some SGE commands and their Slurm equivalents.

| User Command | SGE | Slurm |
|:------------:|-----|-------|
| remote login | `qrsh/qlogin` | `srun --pty bash` |
| run interactively | N/A | `srun --pty program` |
| submit job | `qsub script.sh` | `sbatch script.sh` |
| delete job | `qdel job-id` | `scancel job-id` |
| job status by job id | N/A | `squeue --job job-id` |
| detailed job status | `qstat -u '*' -j job-id` | `sstat job-id` |
| job status of your jobs | `qstat` | `squeue --me` |
| job status by user | `qstat -u user` | `squeue -u user` |
| hold job | `qhold job-id` | `scontrol hold job-id` |
| release job | `qrls job-id` | `scontrol release job-id` |
| queue list | `qconf -sql` | `scontrol show partitions` |
| node list | `qhost` | `sinfo -N` OR `scontrol show nodes` |
| cluster status | `qhost -q` | `sinfo` |
| show node resources | N/A | `sinfo "%n %G"` |


| Job Specification | SGE | Slurm |
|:-----------------:|-----|-------|
| script directive marker | `#$` | `#SBATCH` |
| (run in queue) | `-q queue` | `-p queue` |
| allocated nodes | N/A | `-N min[-max]` |
| allocate cores | `-pe smp count` | `-n count` |
| limit running time | `-l h_rt=time` | `-t days-hh:mm:s` |
| redirectd stdout | `-o file` | `-o file` |
| redirect stderr | `-e file` | `-e file` |
| combine stdout/stderr | `-j yes` | `-o without -e` |
| copy environment | `-V` | `--export=ALL\|NONE\|variables` |
| email notification | `-m abe` | `--mail-type=events` |
| send email to | `-M email` | `--mail-user=email` |
| job name | `-N name` | `--job-name=name` |
| restart job | `-r yes|no` | `--requeue|--no-requeue` |
| working directory | `-wd path` | `--workdir` |
| run exclusively | `-l exclusive` | `--exclusive` OR `--shared` |
| allocate memory | `-l h_vmem=size` | `--mem=mem` OR `--mem-per-cpu=mem` |
| wait for job | `-hold_jid jid` | `--depend state:job` |
| select target host | `-l hostname=host1\|host1` | `--nodelist=nodes` AND/OR `--exclude` |
| allocate GPU | `-l gpu=1` | `--gres=gpu:tesla:count` |
