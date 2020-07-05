# How-To: Run Keras (Multi-GPU)


Because the GPU nodes `med030[1-4]` has four GPU units we can train a model by using multiple GPUs in parallel. This How-To gives an example with Keras 2.2.4 together and tensorflow. Finally soem hints how you can submit a job on the cluster.

!!! hint

    With tensorflow > 2.0 and newer keras version the [`multi_gpu_model`](https://www.tensorflow.org/api_docs/python/tf/keras/utils/multi_gpu_model) is deprecated and you have to use the [`MirroredStrategy`](https://www.tensorflow.org/api_docs/python/tf/distribute/MirroredStrategy).
    
  
# Keras code

we need to import the `multi_gpu_model` model from `keras.utils` and have to pass our actual model (maybe sequential Keras model) into it. In general Keras automatically configures the number of available nodes (`gpus=None`). This seems not to work on our system. So we have to specify the numer of GPUs, e.g. two with `gpus=2`. We put this in a try catch environment that it will also work on CPUs. 

```python
from keras.utils import multi_gpu_model

try: 
    model = multi_gpu_model(model, gpus=2) 
except:
    pass
```

That's it!

Please read [here](../connect/gpu-nodes.md) on how to submit jobs to the GPU nodes.


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
