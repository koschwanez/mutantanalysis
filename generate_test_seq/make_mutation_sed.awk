#! /usr/bin/awk
# parses mutation file to make sed script

BEGIN { FS = "\t" }

$4 ~ r_or_a && $5 ~ seg_perc {
    printf ( "s/%s/%s/g\n", $1, $2 )
}


