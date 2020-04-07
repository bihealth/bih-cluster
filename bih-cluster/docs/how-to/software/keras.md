# How-To: Run Keras (Multi-GPU)

!!! todo "Not yet updated to Slurm"

    TODO: This still needs to be updated to Slurm.

Because the GPU node `med0405` has two GPU units we can train a model by using both GPUs in parallel. This How-To gives an example with Keras 2.2.4 together and tensorflow. Finally soem hints how you can submit a job on the cluster.


# Keras code

we need to import the `multi_gpu_model` model from `keras.utils` and have to pass our actual model (maybe sequential Keras model) into it. In general Keras automatically configures the number of available nodes (`gpus=None`). This seems not to work on our system. So we have to specify them `gpus=2`. We put this in a try catch environment that it will also work on CPUs. 

```python
from keras.utils import multi_gpu_model

try: 
    model = multi_gpu_model(model, gpus=2) 
except:
    pass
```

That's it!

# qsub

If you submit a job to the GPU make sure that you use the `gpu` que and specify both GPUs. We can only allocate one CPU thread in the GPU queue (but more will be used anyway...)

```
-P gpu -pe smp 1 -l gpu=2 h_vmem=100G
```

# Conda environment

All this was tested with the following conda environment:

```yml
name: cuda channels: 
- conda-forge
- bioconda
- defaults
dependencies:
- keras=2.2.4
- python=3.6.7
- tensorboard=1.12.0
- tensorflow=1.12.0
- tensorflow-base=1.12.0
- tensorflow-gpu=1.12.0
```