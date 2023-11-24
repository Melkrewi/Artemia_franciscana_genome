#!/usr/bin/perl
# Written by Will Gammerdinger at IST Austria on September 14th. The original copy was lost in the September 13th incident. Output was md5sum confirmed with the output from the original version. 
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Make variables for GetOptions
my ($blast_output_1, $blast_output_2, $output_file, $help) = ("empty", "empty", "empty");

# Use GetOptions to take entries from the command-line and assign them to variables
GetOptions(
    "blast_output_1=s"          => \$blast_output_1,
    "blast_output_2=s"          => \$blast_output_2,
    "output_file=s"             => \$output_file,
) or Usage ();

# If any of the input variables is undefined then kill the program and trigger this message to show the proper input format and which variables are assigned to currently
if ($blast_output_1 eq "empty" || $blast_output_2 eq "empty" || $output_file eq "empty"){
    die "\nERROR: The format should be perl Reciprocal_BLAST.pl --blast_output_1=blast_output_1_format_6.tab --blast_output_2=blast_output_2_format_6.tab --output_file=output_file.txt\n\nYour current inputs are:\n\nblast_output_1: $blast_output_1\nblast_output_2: $blast_output_2\noutput_file: $output_file\n\n"
}

# Open the first BLAST file
open (my $BLAST_OUTPUT_1, "<$blast_output_1");

# Make hash to hold blast_output_1's BLAST results
my %blast_hash;

# Read each line of blast_output_1
while (my $line_1 = <$BLAST_OUTPUT_1>){
    # Split the line into an array on the tabs
    my @array_of_line_1 = split(/\t/, $line_1);
    # Assign the query to a variable
    my $query_1=$array_of_line_1[0];
    # Assign the hit to a variable
    my $hit_1=$array_of_line_1[1];
    # Read query and hit into the blast_hash
    $blast_hash{$query_1}=$hit_1;
}

# Close the first BLAST file
close $BLAST_OUTPUT_1;

# Open the first BLAST file
open (my $BLAST_OUTPUT_2, "<$blast_output_2");

# Open the first BLAST file
open (my $OUTPUT, ">$output_file");


# Read each line of blast_output_2
while (my $line_2 = <$BLAST_OUTPUT_2>){
    # Split the line into an array on the tabs
    my @array_of_line_2 = split(/\t/, $line_2);
    # Assign the query to a variable
    my $query_2=$array_of_line_2[0];
    # Assign the hit to a variable
    my $hit_2=$array_of_line_2[1];
    # If the hit is a query in the blast hash
    if (exists $blast_hash{$hit_2}){
	# If the hit when used as the query in the blast hash returns the query
	if ($blast_hash{$hit_2} eq $query_2){
	    # Print the accession numbers to the output file
	    print $OUTPUT "$query_2\t$hit_2\n"
	}
    }
}

# Close the second BLAST file
close $BLAST_OUTPUT_2;

# Close the output file
close $OUTPUT;

# If the Usage subroutine is triggered
sub Usage 
{
    # Print an usage example
    printf STDERR "\nThe format should be perl Reciprocal_BLAST.pl --blast_output_1=blast_output_1_format_6.tab --blast_output_2=blast_output_2_format_6.tab --output_file=output_file.txt\n\n";
    # Exit Reciprocal_BLAST.pl
    exit;
}
