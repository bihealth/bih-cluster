# AlphaFold3

the public databases are centrally stored at `/data/cephfs-1/work/projects/alphafold3/public-databases`. 

CUBI does not provide active maintenance, use at your own risk. 
to get access to the project space, write an email to hpc-helpdesk@bih-charite.de.

please apply for your own license [here](https://docs.google.com/forms/d/e/1FAIpQLSfWZAgo1aYk0O4MuAXZj8xRQ8DafeFJnldNOnh_13qAx2ceZw/viewform)  to get the model parameters / weights (~1 GB file)

a singularity image with a working AlphaFold3 installation is available under `apptainer`, 
and there’s a shell script (written mostly by chatGPT) that emulates the `run_alphafold3.py `behavior 
but executes the container under the hood and mounts all directories into it.

so for instance:

```
bash run_alphafold3.sh \
   --json_path path/to/input/input.json \ 
   --model_dir /path/to/alphafold_parameters \ 
   --db_dir /data/cephfs-1/work/projects/alphafold3/public-databases \ 
   --output_dir af_output      
```

note that the V100 compute nodes (`hpc-gpu-[1-8]`) require the extra flag `--flash_attention_implementation=xla`
