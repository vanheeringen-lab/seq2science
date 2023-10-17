# Seq2science profiles

It is possible that you have your own custom settings you always use, e.g. when running seq2science on a cluster.
For these cases we made generic profiles, and can be found here: https://github.com/vanheeringen-lab/seq2science-profiles

Support for profiles is **experimental**, and might not always work as intended. Please let us know if you come across any issues.

Running seq2science after installation of a profile is as easy as:

```
seq2science run {workflow} --profile slurm
```
