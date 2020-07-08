# Convert a time in "days-hh:mm:ss" format or "hh:mm:ss" format into minutes.
function time_in_minutes(time) {
    n = split(time, a, "-");
    if (n > 1) {
        days = a[1];
        time = a[2];
    } else {
        days = 0;
    }

    n = split(time, a, ":");
    if (n == 1) {
        hours = 0;
        minutes = 0;
        seconds = a[1]
    } else if (n == 2) {
        hours = 0;
        minutes = a[1];
        seconds = a[2];
    } else {
        hours = a[1];
        minutes = a[2];
        seconds = a[3];
    }

    return (days * 24 * 60 + hours * 60 + minutes);
}

# Return memory in gigabytes.
function mem_in_gbytes(mem, cpus, mult) {
    if (mult == 0) {
        suffix = substr(mem, length(mem));
        mem = substr(mem, 1, length(mem) - 1);
        mult = (suffix = "n") ? 1 : cpus;
    }

    suffix = substr(mem, length(mem));
    prefix = substr(mem, 1, length(mem) - 1);

    if (suffix == "G") {
        mem = prefix;
    } else if (suffix == "M") {
        mem = prefix / 1024;
    } else if (suffix = "K") {
        mem = prefix / 1024 / 1024;
    }

    return mem * mult;
}

# Initialize.
BEGIN {
    idx_elapsed = 0;
    idx_total = 0;
    idx_cpus = 0;
    idx_reqmem = 0;
    idx_maxrss = 0;
}

# Get field indexes and augment the hader line.
(NR == 1) {
    first_nf = NF;

    for (i = 1; i <= NF; i++) {
        if ($i == "Elapsed") {
            idx_elapsed = i;
        } else if ($i == "TotalCPU") {
            idx_total = i;
        } else if ($i == "AllocCPUS") {
            idx_cpus = i;
        } else if ($i == "ReqMem") {
            idx_reqmem = i;
        } else if ($i == "MaxRSS") {
            idx_maxrss = i;
        }
    }

    if (idx_elapsed == 0) {
        print "Header \"Elapsed\" not found!" >"/dev/stderr"
        exit 1;
    }
    if (idx_total == 0) {
        print "Header \"TotalCPU\" not found!" >"/dev/stderr"
        exit 1;
    }
    if (idx_cpus == 0) {
        print "Header \"AllocCPUS\" not found!" >"/dev/stderr"
        exit 1;
    }
    if (idx_reqmem == 0) {
        print "Header \"ReqMem\" not found!" >"/dev/stderr"
        exit 1;
    }
    if (idx_maxrss == 0) {
        print "Header \"MaxRSS\" not found!" >"/dev/stderr"
        exit 1;
    }

    printf("%s%s\n", $0, "    EmpPar   ParEff  DiffMEM");
}
# Print the augmented header line.
(NR == 2) {
    printf("%s%s\n", $0, " --------- -------- --------");
}
# Augment the row
(NR > 2) {
    elapsed = time_in_minutes($idx_elapsed);
    total = time_in_minutes($idx_total);
    cpus = $idx_cpus;
    reqmem = mem_in_gbytes($idx_reqmem, cpus, 0);
    maxrss = (first_nf != NF) ? 0 : mem_in_gbytes($idx_maxrss, cpus, 1);
    emp_par = (total == 0) ? 0 : total / elapsed;
    eff_par = (total == 0) ? 0 : total / (elapsed * cpus);
    if (maxrss == 0) {
        printf("%s %8.2f %8.2f %8s\n", $0, emp_par, eff_par, "-");
    } else {
        printf("%s %8.2f %8.2f %8.2f\n", $0, emp_par, eff_par, reqmem - maxrss);
    }
}