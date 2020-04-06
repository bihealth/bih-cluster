# How-To: Setup TensorFlow

TensorFlow is a package for deep learning with optional support for GPUs.
You can find the original TensorFlow installation instructions [here](https://www.tensorflow.org/versions/r1.3/get_started/os_setup.html#download-and-setup).

This article describes how to set up TensorFlow with GPU support using Conda.
This how-to assumes that you have just connected to a GPU node via `qrsh -P gpu -l gpu=1` (with `-P gpu`, the `gpu` SGE project gives you access to the GPU node and `-l gpu=1` reserves the GPU resource) .

At the time of writing, Tensorflow was available with GPU support from conda in version 1.2.

:lower_left_ballpoint_pen: You're welcome to update this (and any other) wiki page with new information.
This is the only way to keep the Wiki up to date.

This tutorial assumes, that conda has been set up as described in [Software Management](Manual-Software-Management/home).

## Create conda environment

```terminal
$ conda create -n cuda python=3.6
$ source activate cuda
$ conda install tensorflow-gpu==1.12.0
```

## Run TensorFlow example

We now trying to run the example from the TensorFlow documentation.

:warning: Note that the two `export` lines are important. :exclamation:

```terminal
## these two lines are important
$ export CUDA_HOME=$(dirname $(dirname $(which python)))
$ export LD_LIBRARY_PATH=$CUDA_HOME/lib:$LD_LIBRARY_PATH

$ python
>>> import tensorflow as tf
>>> hello = tf.constant('Hello, TensorFlow!')
>>> sess = tf.Session()
2017-09-19 18:07:40.571996: W tensorflow/core/platform/cpu_feature_guard.cc:45] The TensorFlow library wasn't compiled to use SSE4.1 instructions, but these are available on your machine and could speed up CPU computations.
2017-09-19 18:07:40.572090: W tensorflow/core/platform/cpu_feature_guard.cc:45] The TensorFlow library wasn't compiled to use SSE4.2 instructions, but these are available on your machine and could speed up CPU computations.
2017-09-19 18:07:40.572129: W tensorflow/core/platform/cpu_feature_guard.cc:45] The TensorFlow library wasn't compiled to use AVX instructions, but these are available on your machine and could speed up CPU computations.
2017-09-19 18:07:43.742453: I tensorflow/core/common_runtime/gpu/gpu_device.cc:940] Found device 0 with properties:
name: Tesla K20Xm
major: 3 minor: 5 memoryClockRate (GHz) 0.732
pciBusID 0000:04:00.0
Total memory: 5.57GiB
Free memory: 5.50GiB
2017-09-19 18:07:43.983863: W tensorflow/stream_executor/cuda/cuda_driver.cc:523] A non-primary context 0x3b11f70 exists before initializing the StreamExecutor. We haven't verified StreamExecutor works with that.
2017-09-19 18:07:43.985238: I tensorflow/core/common_runtime/gpu/gpu_device.cc:940] Found device 1 with properties:
name: Tesla K20Xm
major: 3 minor: 5 memoryClockRate (GHz) 0.732
pciBusID 0000:42:00.0
Total memory: 5.57GiB
Free memory: 5.50GiB
2017-09-19 18:07:43.985403: I tensorflow/core/common_runtime/gpu/gpu_device.cc:832] Peer access not supported between device ordinals 0 and 1
2017-09-19 18:07:43.985440: I tensorflow/core/common_runtime/gpu/gpu_device.cc:832] Peer access not supported between device ordinals 1 and 0
2017-09-19 18:07:43.985498: I tensorflow/core/common_runtime/gpu/gpu_device.cc:961] DMA: 0 1
2017-09-19 18:07:43.985517: I tensorflow/core/common_runtime/gpu/gpu_device.cc:971] 0:   Y N
2017-09-19 18:07:43.985533: I tensorflow/core/common_runtime/gpu/gpu_device.cc:971] 1:   N Y
2017-09-19 18:07:43.985657: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1030] Creating TensorFlow device (/gpu:0) -> (device: 0, name: Tesla K20Xm, pci bus id: 0000:04:00.0)
2017-09-19 18:07:43.985683: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1030] Creating TensorFlow device (/gpu:1) -> (device: 1, name: Tesla K20Xm, pci bus id: 0000:42:00.0)
>>> print(sess.run(hello))
b'Hello, TensorFlow!'
>>> a = tf.constant(10)
>>> b = tf.constant(32)
>>> print(sess.run(a + b))
42
>>> exit()
```

Check that you get the lines metioning `Tesla K20Xm` and the two lines mentioning `gpu:0` and `gpu:1`.
