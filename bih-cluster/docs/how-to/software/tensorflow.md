# How-To: Setup TensorFlow

TensorFlow is a package for deep learning with optional support for GPUs.
You can find the original TensorFlow installation instructions [here](https://www.tensorflow.org/install).

This article describes how to set up TensorFlow with GPU support using Conda.
This how-to assumes that you have just connected to a GPU node via `srun --mem=10g --partition=gpu --gres=gpu:tesla:1 --pty bash -i`.
Note that you will need to allocate "enough" memory, otherwise your python session will be `Killed` because of too little memory.
You should read the [How-To: Connect to GPU Nodes](../../how-to/connect/gpu-nodes/) tutorial on an explanation of how to do this and to learn how to register for GPU usage.

This tutorial assumes, that conda has been set up as described in [Software Management]((../../best-practice/software-installation-with-conda.md).

## Create conda environment

We recommend that you install mamba first with `conda install -y mamba` and use this C++ reimplementation of the `conda command` as follows.

```terminal
$ conda create -y -n python-tf tensorflow-gpu
$ conda activate python-tf
```

Let us verify that we have Python and TensorFlow installed.
You might get different versions you could pin the version on installing with `conda create -y -n python-tf python==3.9.10 tensorflow-gpu==2.6.2

```terminal
$ python --version
Python 3.9.10
$ python -c 'import tensorflow; print(tensorflow.__version__)'
2.6.2
```

We thus end up with an installation of Python 3.9.10 with tensorflow 2.6.2.

## Run TensorFlow Example

Let us now see whether TensorFlow has recognized our GPU correctly.

```terminal
$ python
>>> import tensorflow as tf
>>> print("TensorFlow version:", tf.__version__)
TensorFlow version: 2.6.2
>>> print(tf.config.list_physical_devices())
[PhysicalDevice(name='/physical_device:CPU:0', device_type='CPU'), PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')]
```

Yay, we can proceed to run the [Quickstart Tutorial](https://www.tensorflow.org/tutorials/quickstart/beginner).

```
>>> mnist = tf.keras.datasets.mnist
>>> (x_train, y_train), (x_test, y_test) = mnist.load_data()
>>> x_train, x_test = x_train / 255.0, x_test / 255.0
>>> model = tf.keras.models.Sequential([
...   tf.keras.layers.Flatten(input_shape=(28, 28)),
...   tf.keras.layers.Dense(128, activation='relu'),
...   tf.keras.layers.Dropout(0.2),
...   tf.keras.layers.Dense(10)
... ])
>>> predictions = model(x_train[:1]).numpy()
>>> predictions
array([[-0.50569224,  0.26386747,  0.43226188,  0.61226094,  0.09630793,
         0.34400576,  0.9819117 , -0.3693726 ,  0.5221357 ,  0.3323232 ]],
      dtype=float32)
>>> tf.nn.softmax(predictions).numpy()
array([[0.04234391, 0.09141268, 0.10817807, 0.12951255, 0.07731011,
        0.09903987, 0.18743432, 0.04852816, 0.11835073, 0.09788957]],
      dtype=float32)
>>> loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
>>> loss_fn(y_train[:1], predictions).numpy()
2.3122327
>>> model.compile(optimizer='adam',
...               loss=loss_fn,
...               metrics=['accuracy'])
>>> model.fit(x_train, y_train, epochs=5)
2022-03-09 17:53:47.237997: I tensorflow/compiler/mlir/mlir_graph_optimization_pass.cc:185] None of the MLIR Optimization Passes are enabled (registered 2)
Epoch 1/5
1875/1875 [==============================] - 3s 1ms/step - loss: 0.2918 - accuracy: 0.9151
Epoch 2/5
1875/1875 [==============================] - 3s 1ms/step - loss: 0.1444 - accuracy: 0.9561
Epoch 3/5
1875/1875 [==============================] - 3s 1ms/step - loss: 0.1082 - accuracy: 0.9674
Epoch 4/5
1875/1875 [==============================] - 3s 1ms/step - loss: 0.0898 - accuracy: 0.9720
Epoch 5/5
1875/1875 [==============================] - 3s 1ms/step - loss: 0.0773 - accuracy: 0.9756
<keras.callbacks.History object at 0x154e81360190>
>>> model.evaluate(x_test,  y_test, verbose=2)
313/313 - 0s - loss: 0.0713 - accuracy: 0.9785
[0.0713074803352356, 0.9785000085830688]
>>> probability_model = tf.keras.Sequential([
...   model,
...   tf.keras.layers.Softmax()
... ])
>>> probability_model(x_test[:5])
<tf.Tensor: shape=(5, 10), dtype=float32, numpy=
array([[1.2339272e-06, 6.5599060e-10, 1.0560590e-06, 5.9356184e-06,
        5.3691075e-12, 1.4447859e-07, 5.4218874e-13, 9.9996936e-01,
        1.0347234e-07, 2.2147648e-05],
       [2.9887938e-06, 6.8461006e-05, 9.9991941e-01, 7.2003731e-06,
        2.9751782e-13, 8.2818183e-08, 1.4307782e-06, 2.3203837e-13,
        4.7433215e-07, 2.9504194e-14],
       [1.8058477e-06, 9.9928612e-01, 7.8716243e-05, 3.9140195e-06,
        3.0842333e-05, 9.4537208e-06, 2.2774333e-05, 4.5549971e-04,
        1.1015874e-04, 6.9138093e-07],
       [9.9978787e-01, 3.0206781e-08, 2.8528208e-05, 8.5581682e-08,
        1.3851340e-07, 2.3634559e-06, 1.8480707e-05, 1.0153375e-04,
        1.1583331e-07, 6.0887167e-05],
       [6.4914235e-07, 2.5808356e-08, 1.8225538e-06, 2.3215563e-09,
        9.9588013e-01, 4.6049720e-08, 3.8903639e-07, 2.9772724e-05,
        4.3141077e-07, 4.0867776e-03]], dtype=float32)>
>>> exit()
```

## Writing TensorFlow Slurm Jobs

Writing Slurm jobs using TensorFlow is as easy as creating the following scripts.

`tf_script.py`

```python
#/usr/bin/env python

import tensorflow as tf
print("TensorFlow version:", tf.__version__)
print(tf.config.list_physical_devices())

mnist = tf.keras.datasets.mnist

(x_train, y_train), (x_test, y_test) = mnist.load_data()
x_train, x_test = x_train / 255.0, x_test / 255.0


model = tf.keras.models.Sequential([
  tf.keras.layers.Flatten(input_shape=(28, 28)),
  tf.keras.layers.Dense(128, activation='relu'),
  tf.keras.layers.Dropout(0.2),
  tf.keras.layers.Dense(10)
])

predictions = model(x_train[:1]).numpy()
print(predictions)

print(tf.nn.softmax(predictions).numpy())

# ... and so on ;-)
```

`tf_job.sh`

```bash
#!/usr/bin/bash

#SBATCH --job-name=tf-job
#SBATCH --mem=10g
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tesla:1

source $HOME/work/miniconda3/bin/activate
conda activate python-tf

python tf_script.py &>tf-out.txt
```

And then calling

```terminal
$ sbatch tf_job.sh
```

You can find the reuslts in `tf-out.txt` after completion.

```
$ cat tf-out.txt 
2022-03-09 18:05:54.628846: I tensorflow/core/platform/cpu_feature_guard.cc:142] This TensorFlow binary is optimized with oneAPI Deep Neural Network Library (oneDNN) to use the following CPU instructions in performance-critical operations:  SSE4.1 SSE4.2 AVX AVX2 AVX512F FMA
To enable them in other operations, rebuild TensorFlow with the appropriate compiler flags.
2022-03-09 18:05:56.999848: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1510] Created device /job:localhost/replica:0/task:0/device:GPU:0 with 30988 MB memory:  -> device: 0, name: Tesla V100-SXM2-32GB, pci bus id: 0000:18:00.0, compute capability: 7.0
TensorFlow version: 2.6.2
[PhysicalDevice(name='/physical_device:CPU:0', device_type='CPU'), PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')]
[[-0.07757086  0.04676083  0.9420195  -0.59902835 -0.26286742 -0.392514
   0.3231195  -0.17169198  0.3480805   0.37013203]]
[[0.07963609 0.09017922 0.22075593 0.04727634 0.06616627 0.05812084
  0.11888511 0.07248258 0.12188996 0.12460768]]
```
