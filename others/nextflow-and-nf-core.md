---
description: >-
  Nextflow is an incredibly powerful and flexible workflow language, and nf-core
  is a curated set of open‑source analysis pipelines built using Nextflow.
---

# nextflow & nf-core

## Links

* Nextflow: [https://www.nextflow.io/](https://www.nextflow.io/)
  * Tutorial: [https://training.nextflow.io/basic\_training/](https://training.nextflow.io/basic_training/)
* nf-core: [https://nf-co.re/](https://nf-co.re/)

## Operating nf-core offline in hpc

reference: [https://nf-co.re/docs/usage/offline](https://nf-co.re/docs/usage/offline)

```sh
nf-core download nf-core/atacseq --revision 2.1.2 --container-system singularity

tar -xvzf <file.tar.gz> -C <some_custom_folder_name>
```

* unzipped folder should contain:
  * workflow (the pipeline files)
  * config (a copy of nf-core/configs)
  * singularity (if —container singularity is used)
*   to run the pipeline

    * <mark style="color:red;">`nextflow run <download_directory>/workflow [pipeline flags]`</mark>

    ```bash
    ./nextflow run ./nf-core-cutandrun_3.2.1/3_2_1/ -profile "test"
    ```
*   had to install the newer version of nextflow to install nf-validation <mark style="color:red;">`./nextflow plugin install nf-validation`</mark>

    ```bash
    module load jre/jdk-13

    curl -s <https://get.nextflow.io> | bash

    ## by using the executable file of nextflow directly instead of the installed version
    ./nextflow plugin install nf-validation
    ```
*   setting the cache directory

    ```bash
    export NXF_SINGULARITY_CACHEDIR=/hpctmp/userid/nf-core-cutandrun_3.2.1/singularity-images/
    ```
* nextflow profile configuration: [https://github.com/nf-core/configs](https://github.com/nf-core/configs)
  * [https://www.nextflow.io/docs/latest/config.html#scope-singularity](https://www.nextflow.io/docs/latest/config.html#scope-singularity)

{% code title="hpc.config" %}
```
## example of hpc configuration
params {
   config_profile_description = 'NUS HPC'
   config_profile_contact = 'User'

   max_memory = 128.GB
   max_cpus = 20
   max_time = 1000.h
}

singularity {
   enabled = true
   cacheDir = "/scratch2/userid/nf-core-atacseq_2.1.2/singularity-images/"
   autoMounts = true
}

plugins {
  id 'nf-validation@1.1.3'
}
```
{% endcode %}

## nextflow - tracing

To show information about executed pipelines in the current folder: [https://www.nextflow.io/docs/latest/tracing.html](https://www.nextflow.io/docs/latest/tracing.html)

```
nextflow log <run name> [options]
```
