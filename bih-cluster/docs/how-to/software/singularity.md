# Using Singularity (with Docker Images)

Singularity (https://sylabs.io/docs/) is a popular alternative to docker, because it does not require to run as a privileged user.
Singualrity can run Docker images out-of-the-box by converting them to the singularity image format.
The following guide gives a quick dive into using docker images with singularity.

## Quickstart

!!! important "Link ~/.singularity to ~/work/.singularity"

    Because you only have a quota of 1 GB in your home directory, you should symlink `~/.singularity` to `~/work/.singularity`.

    ```bash
    host:~$ mkdir -p ~/work/.singularity && ln -sr ~/work/.singularity ~/.singularity
    ```

    In case you already have a singularity directory:

    ```bash
    host:~$ mv ~/.singularity ~/work/.singularity && ln -sr ~/work/.singularity ~/.singularity
    ```

Run a bash in a docker image:

```bash
host:~$ singularity bash docker://godlovedc/lolcow
```

Run a command in a docker image:

```bash
host:~$ singularity exec docker://godlovedc/lolcow echo "hello, hello!"
```

Run a bash in a docker image, enable access to the cuda driver (--nv) and mount a path (--bind or -B):

```bash
host:~$ singularity bash --nv --bind /path_on_host/:/path_inside_container/ docker://godlovedc/lolcow
```

## Some Caveats and Notes

**Caveats**

- The default singularity images format (.sif) is read-only.
- By default singularity mounts /home/$USER, /tmp, and $PWD in the container.

**Notes**

- Environment variables can be provided by setting them in the bash and adding the prefix `SINGULARITYENV_`:
    ```bash
    host:~$ SINGULARITYENV_HELLO=123 singularity bash docker://godlovedc/lolcow echo $HELLO
    ```
- Calling `singularity shell` or `singularity exec` uses as cwd the host callers cwd not the one set in the Dockerfile.
  One can change this by setting `--pwd`.

## Referencing/Providing Docker Images

### Option 1: Use Docker Images via Docker Hub

The easiest variant to run a docker image available via a docker hub is by specifying its url.
This causes singularity to download the image and convert it to a singularity image:

```bash
host:~$ singularity run docker://godlovedc/lolcow
```

or to open a shell inside the image

```bash
host:~$ singularity bash docker://godlovedc/lolcow
```

Furthermore, similar to docker, one can pull (and convert) remote image with the following call:

```bash
host:~$ singularity pull docker://godlovedc/lolcow
```

In case your registry requires authentication you can provide it via a prompt by adding the option `--docker-login`:

```bash
host:~$ singularity pull --docker-login docker://ilumb/mylolcow
```

or by setting the following environment variables:

```bash
host:~$ export SINGULARITY_DOCKER_USERNAME=ilumb
host:~$ export SINGULARITY_DOCKER_PASSWORD=<redacted>
host:~$ singularity pull docker://ilumb/mylolcow
```

More details can be found [in the Singularity documentation](https://sylabs.io/guides/3.5/user-guide/singularity_and_docker.html).

### Option 2: Converting Docker Images

Another option is to convert your docker image into the Singularity image format.
This can be easily done using the docker images provided by [docker2singularity](https://github.com/singularityhub/docker2singularity).

To convert the docker image `docker_image_name` to the singularity image `singularity_image_name` one can use the following command line.
The output image will be located in `output_directory_for_images`.

```bash
host:~$ docker run -v /var/run/docker.sock:/var/run/docker.sock --privileged -t --rm quay.io/singularity/docker2singularity -v /output_directory_for_images/:/output --name singularity_image_name docker_image_name
```

The resulting image can then directly be used as image:

```bash
host:~$ singularity bash singularity_image_name
```

## Conversion Compatibility

Here are some tips for making Docker images compatible with Singularity taken from [docker2singulrity](https://github.com/singularityhub/docker2singularity):

- Define all environmental variables using the ENV instruction set. Do not rely on `~/.bashrc`, `~/.profile`, etc.
- Define an `ENTRYPOINT` instruction set pointing to the command line interface to your pipeline.
- Do not define `CMD` - rely only on `ENTRYPOINT`.
- You can interactively test the software inside the container by overriding the `ENTRYPOINT docker run -i -t --entrypoint /bin/bash bids/example`.
- Do not rely on being able to write anywhere other than the home folder and /scratch.
  Make sure your container runs with the `--read-only --tmpfs /run --tmpfs /tmp parameters` (this emulates the read-only behavior of Singularity).
- Don't rely on having elevated user permissions.
- Don't use the `USER` instruction set.
