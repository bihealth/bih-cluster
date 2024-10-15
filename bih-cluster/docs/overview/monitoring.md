# Monitoring

We currently provide you only with Ganglia for monitoring the cluster status.

!!! cite "Ganglia Ganglia was shut down in favour of HPC Metrics"
    [Ganglia](https://metrics.cubi.bihealth.org/public-dashboards/dc3e4d5b1ea049429abf39e412c47302?orgId=1&refresh=1m)

## Using Ganglia

Go to the following address and login with your home organization (Charite or MDC):

- https://hpc-ganglia.cubi.bihealth.org


!!! cite "Ganglia does not know about Slurm"
    Ganglia will not show you anything about the Slurm job schedulign system.
    If a job uses a whole node but uses no CPUs then this will be displayed as unused in Ganglia.
    However, Slurm would not schedule another job on this node.

You will be show a screen as shown below.
This allows you to get a good idea of what is going on on the HPC.

![](figures/Ganglia_Example.png)

By default you will be shown the cluster usage of the last day.
You can quickly switch to report for two or four hours as well, etc.

In the first row of pictures you see the number of total CPUs (actually hardware threads), number of hosts seen as up and down by Ganglia, and cluster load/utilization.
You will then see the overall cluster load, memory usage, CPU usage, and network utilization across the selected time period.

!!! cite "Linux load is not intuitive"
    Note that the technical details behind Linux **load** is not very interactive.
    It is incorporating much more than just the CPU usage.
    You can find a quite comprehensive [treatement of Linux Load here](https://www.brendangregg.com/blog/2017-08-08/linux-load-averages.html).

We are using a fast shared storage system and almost no local storage (except in `/tmp`).
Also, almost no jobs use MPI or other heavy network communication.
Thus, the network utilization is a good measure of the I/O on the cluster.

Below, you can drill down into various metrics and visualize them historically.
Just try it out and find your way around, you cannot break anything.
Sadly, there is no good documentation of Ganglia online.

## Aggregate GPU Utilization Visualization

Ganglia allows you to obtain metrics in several interesting and useful ways.
If you click on "Aggregate Graphs" then you could enter the following values to get an overview of the live GPU utilization.

- Title: `Aggreate GPU Utilization`
- Host Regular expression: `hpc-gpu-.*`
- Metric Regular Expressions: `gpu._util`
- Graph Type: `Stacked`
- Legend Options: `Hide legend`

Then click `Create Graph`.

![](figures/Ganglia_Aggregate_GPUs.png)

If a GPU is fully used, it will contribute 100 points on the vertical axis.
See above for an example, and here is a direct link:

- [Aggregate GPU Utilization](https://hpc-ganglia.cubi.bihealth.org/ganglia/graph_all_periods.php?title=Aggregate+GPU+Utilization&cs=&ce=&vl=&x=&n=&hreg%5B%5D=hpc-gpu-.*&mreg%5B%5D=gpu._util&gtype=stack&glegend=hide&aggregate=1)
