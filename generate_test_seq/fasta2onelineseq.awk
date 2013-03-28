#! /usr/bin/awk
# Inspired by blossomassociates.net fasta2seq.awk script

#description line
/^>/ { 
    if (current_sequence) {
        print current_sequence
        current_sequence = ""
    }
    print $0 
}

# comment line
/^;/ {print $0}

# collect sequence
/^[^>;]/{ current_sequence = current_sequence $0 }

END { if (current_sequence) print current_sequence }
