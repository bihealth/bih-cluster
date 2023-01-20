# Using Apptainer (with Docker Images)

!!! note

    Singularity is now Apptainer!
    While Apptainer provides an `singularity` alias for backwards compatibility,
    it is recommanded to adapt all workflows to use the new binary `apptainer`.

Apptainer (https://apptainer.org/) is a popular alternative to docker, because it does not require to run as a privileged user.
Apptainer can run Docker images out-of-the-box by converting them to the apptainer image format.
The following guide gives a quick dive into using docker images with apptainer.

!!! important "Build on your workstation, run on the HPC"

    Building images using Apptainer requires root privileges.
    We cannot give you these permissions on the BIH HPC.
    Thus, you will have to build the images on your local workstation (or anywhere where you have root access).
    You can then run the built images on the BIH HPC.

    This is also true for the `--writeable` flag.
    Apparently it needs root permissions which you don't have on the cluster.

## Quickstart

!!! important "Link ~/.apptainer to ~/work/.apptainer"

    Because you only have a quota of 1 GB in your home directory, you should symlink `~/.apptainer` to `~/work/.apptainer`.

    ```bash
    host:~$ mkdir -p ~/work/.apptainer && ln -sr ~/work/.apptainer ~/.apptainer
    ```

    In case you already have a apptainer directory:

    ```bash
    host:~$ mv ~/.apptainer ~/work/.apptainer && ln -sr ~/work/.apptainer ~/.apptainer
    ```

Run a bash in a docker image:

```bash
host:~$ apptainer bash docker://godlovedc/lolcow
```

Run a command in a docker image:

```bash
host:~$ apptainer exec docker://godlovedc/lolcow echo "hello, hello!"
```

Run a bash in a docker image, enable access to the cuda driver (--nv) and mount a path (--bind or -B):

```bash
host:~$ apptainer bash --nv --bind /path_on_host/:/path_inside_container/ docker://godlovedc/lolcow
```

## Some Caveats and Notes

**Caveats**

- The default apptainer images format (.sif) is read-only.
- By default apptainer mounts /home/$USER, /tmp, and $PWD in the container.

**Notes**

- Environment variables can be provided by setting them in the bash and adding the prefix `APPTAINERENV_`:
    ```bash
    host:~$ APPTAINERENV_HELLO=123 apptainer bash docker://godlovedc/lolcow echo $HELLO
    ```
- Calling `apptainer shell` or `apptainer exec` uses as cwd the host callers cwd not the one set in the Dockerfile.
  One can change this by setting `--pwd`.

## Referencing/Providing Docker Images

### Option 1: Use Docker Images via Docker Hub

The easiest variant to run a docker image available via a docker hub is by specifying its url.
This causes apptainer to download the image and convert it to a apptainer image:

```bash
host:~$ apptainer run docker://godlovedc/lolcow
```

or to open a shell inside the image

```bash
host:~$ apptainer bash docker://godlovedc/lolcow
```

Furthermore, similar to docker, one can pull (and convert) remote image with the following call:

```bash
host:~$ apptainer pull docker://godlovedc/lolcow
```

In case your registry requires authentication you can provide it via a prompt by adding the option `--docker-login`:

```bash
host:~$ apptainer pull --docker-login docker://ilumb/mylolcow
```

or by setting the following environment variables:

```bash
host:~$ export APPTAINER_DOCKER_USERNAME=ilumb
host:~$ export APPTAINER_DOCKER_PASSWORD=<redacted>
host:~$ apptainer pull docker://ilumb/mylolcow
```

More details can be found [in the Apptainer documentation](https://apptainer.org/docs/user/main/docker_and_oci.html).

### Option 2: Converting Docker Images

Another option is to convert your docker image into the Apptainer/Singularity image format.
This can be easily done using the docker images provided by [docker2singularity](https://github.com/singularityhub/docker2singularity).

To convert the docker image `docker_image_name` to the apptainer image `apptainer_image_name` one can use the following command line.
The output image will be located in `output_directory_for_images`.

```bash
host:~$ docker run -v /var/run/docker.sock:/var/run/docker.sock --privileged -t --rm quay.io/singularity/docker2singularity -v /output_directory_for_images/:/output --name apptainer_image_name docker_image_name
```

The resulting image can then directly be used as image:

```bash
host:~$ apptainer bash apptainer_image_name
```

## Conversion Compatibility

Here are some tips for making Docker images compatible with Apptainer taken from [docker2singulrity](https://github.com/singularityhub/docker2singularity):

- Define all environmental variables using the ENV instruction set. Do not rely on `~/.bashrc`, `~/.profile`, etc.
- Define an `ENTRYPOINT` instruction set pointing to the command line interface to your pipeline.
- Do not define `CMD` - rely only on `ENTRYPOINT`.
- You can interactively test the software inside the container by overriding the `ENTRYPOINT docker run -i -t --entrypoint /bin/bash bids/example`.
- Do not rely on being able to write anywhere other than the home folder and /scratch.
  Make sure your container runs with the `--read-only --tmpfs /run --tmpfs /tmp parameters` (this emulates the read-only behavior of Apptainer).
- Don't rely on having elevated user permissions.
- Don't use the `USER` instruction set.
