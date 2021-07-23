# Perl-tools

Perl scripts used for different applications.

# GC_percent_window.pl
It calculates the GC content of a multifasta file using a sliding window of predefined length, and printing the results to STDOUT.

  Usage:
  GC_percent_window.pl -i <fasta file> -w <window_size [120bp]>

# bin_depth_values.pl
It bins read depth values from a samtools depth file for plotting.
  
  Usage:
  bin_depth_values.pl -i <depth_file> -b <bin size in bp [120]> -r <rename toxo ctg IDs T/F [F]>

# binner.pl
It groups numerical values into bins specified by the user.

# match_mismatch.pl
It identifies key values that are present or absent from a user-provided list of keys.

  Usage:
  match_mismatch.pl -k <keys_file> [-ki <key file index [0]> ; -d <data_file> | STDIN ; -di <data index [0]> ;-m M|MM [match|mismatch [MM]] ; -h help]

# nr.pl
It removes redundant lines from a file.
  
  Usage:
  nr.pl < redundant_file > non-redundant_file

# stats_her.pl
It calculate several statistics (e.g. mean, median, SD, VAR, sum, max, min, counts) from a numerical vector that is entered from the STDIN

# get_sequence_from_ncbi_v2.sh
It retrieve sequences in bulk from NCBI using a list of accession or gi numbers.
