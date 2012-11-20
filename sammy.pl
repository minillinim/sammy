#!/usr/bin/env perl
###############################################################################
#
#    sammy.pl
#    
#    Stupidly simple paired read simulator - NO ERRORS!
#
#    Copyright (C) Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;

#CPAN modules
use Bio::SeqIO;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands 
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
#printAtStart();
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################
my $global_insert = overrideDefault(500, 'insert');
my $global_stdev = overrideDefault(50, 'deviation');
my $global_orientation = overrideDefault('in', 'orientation');
my $global_read_length = overrideDefault(100, 'length');
my $global_ref_cutoff = overrideDefault(1000, 'cutoff');
my %global_seq_hash = ();       # seqs to cut reads from
my %global_length_hash = ();    # lengths of said seqs
my %global_numbers_hash = ();   # how many reads we'll cut from each

# get the insert sizes
my $global_ins_upper_rand = 0;
my $global_prob_resolution = 10000;
my %global_ins_values = prepInsert();

# open fasta file for all the contigs for reading
my $total_length = 0;
my $seqio = Bio::SeqIO->new(-file => $global_options->{'reference'}, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
    my $header = $seq->id;
    my $sequence = $seq->seq;
    my $seq_length = length($sequence);
    if($seq_length >= $global_ref_cutoff)
    { 
        $global_seq_hash{$header} = $sequence;
        $global_length_hash{$header} = $seq_length;
        $total_length += $seq_length;
    }
}

# work out how many reads to cut from each guy
for my $id (keys(%global_length_hash))
{
    $global_numbers_hash{$id} = int($global_length_hash{$id}/$total_length * $global_options->{'num_reads'});
}

# now make the simulated reads
my $global_seq_length = 0;
my $global_sequence = "";
my $global_rev_sequence = "";
my $read_counter = 1;
foreach my $id (keys(%global_seq_hash))
{
    my $num_reads = $global_numbers_hash{$id};
    
    # set these globals now
    $global_seq_length = $global_length_hash{$id};
    $global_sequence = $global_seq_hash{$id};
    $global_rev_sequence = revcompl($global_sequence);
    
    while($num_reads >= 0){
        my $insert = getInsert($global_ins_upper_rand, $global_prob_resolution);
        my $strand = int(rand(1.99999));
        my ($read_1, $read_2) = cutReads($insert, $strand);
        print ">$read_counter"."_1\n$read_1\n";
        print ">$read_counter"."_2\n$read_2\n";
        $num_reads--;
        $read_counter++;
    } 
}

######################################################################
# CUSTOM SUBS
######################################################################

sub revcompl {
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return scalar reverse $seq;
}

sub cutReads {
    #-----
    # actually cut the reads
    #
    my($insert, $strand) = @_;
    my $read_1 = "";
    my $read_2 = "";

    my $start = int(rand($global_seq_length - $insert));
    
    if($global_orientation eq 'out') {
        # <-- -->
        if($strand == 0) { $read_1 = revcompl(substr($global_sequence, $start, $global_read_length));  $read_2 = substr($global_sequence, $start + $insert - $global_read_length, $global_read_length); }
        else { $read_1 = revcompl(substr($global_rev_sequence, $start, $global_read_length));  $read_2 = substr($global_rev_sequence, $start + $insert - $global_read_length, $global_read_length); }
    }
    elsif($global_orientation eq 'same') {
        # <-- <-- or --> -->
        if($strand == 0) { $read_1 = revcompl(substr($global_sequence, $start, $global_read_length));  $read_2 = revcompl(substr($global_sequence, $start + $insert - $global_read_length, $global_read_length)); }
        else { $read_1 = revcompl(substr($global_rev_sequence, $start, $global_read_length));  $read_2 = revcompl(substr($global_rev_sequence, $start + $insert - $global_read_length, $global_read_length)); }
    }
    elsif($global_orientation eq 'in') {
        # --> <--
        if($strand == 0) { $read_1 = substr($global_sequence, $start, $global_read_length);  $read_2 = revcompl(substr($global_sequence, $start + $insert - $global_read_length, $global_read_length)); }
        else { $read_1 = substr($global_rev_sequence, $start, $global_read_length);  $read_2 = revcompl(substr($global_rev_sequence, $start + $insert - $global_read_length, $global_read_length)); }
    }
    
    return ($read_1, $read_2);
}

sub getInsert {
    #-----
    # Returns an random integer between [ mean insert - (3 * stdev) --> mean insert + (3 * stdev) ]    
    # if run repeatedly, the integers returned will fit a normal distribution
    # this is not 100% true of short read data but it will do    
    # if you're not happy, write your own version here...
    my ($global_ins_upper_rand, $global_prob_resolution) = @_; 
    my $rand_num = 1;
    while($rand_num > $global_ins_upper_rand)
    {  
        $rand_num = int(rand() * $global_prob_resolution) / $global_prob_resolution;
    } 
    return $global_ins_values{$rand_num};
}

sub prepInsert {
    #-----
    # prep the arrays we'll need to use in get insert sizes
    # I wrote this in 2009, I forgot how it works!
    # 
    # Returns a hash of type: float -> int
    #
    # Where float is a number between 0->1
    # 

    my %ins_values = ();
    my $inv_prob_res = 0.0001;


    # these are the limits for the insert sizes
    my $upper_insert_limit = $global_insert + 3 * $global_stdev;
    my $lower_insert_limit = $global_insert - 3 * $global_stdev;


    my $pi = 4 * atan2(1, 1);
    my $cumulative_area = 0;
    my $multiplier = 1 / ($global_stdev * sqrt(2 * $pi));
    my $two_sig_squared =  $global_stdev * $global_stdev * 2;
    my $step = 0.001;
    
    my $x  = 0;
    my $key = $lower_insert_limit; 
    my $prev_cum_area = 0;
    
    # do the first half notch 
    for($x = $lower_insert_limit; $x < $lower_insert_limit + 0.5; $x = $x + $step)
    {
        my $height = $multiplier * exp(-1 * ($x - $global_insert) * ($x - $global_insert) / $two_sig_squared);
        $cumulative_area = $cumulative_area + $height * $step; 
    }
    # update the hash
    my $new_cum_area = int($cumulative_area * $global_prob_resolution);
    for(my $counter = $new_cum_area; $counter >= 0; $counter--)
    {
        $ins_values{($counter / $global_prob_resolution)} = $key;
    }
    $prev_cum_area = $new_cum_area;
    $key++; 
    
    # do all the middle ones
    $lower_insert_limit = $lower_insert_limit + 0.5;
    while($x < $upper_insert_limit - 0.5)
    {  
        for($x = $lower_insert_limit; $x < $lower_insert_limit + 1; $x = $x + $step)
        {
            my $height = $multiplier * exp(-1 * ($x - $global_insert) * ($x - $global_insert) / $two_sig_squared);
            $cumulative_area = $cumulative_area + $height * $step; 
        }
        $lower_insert_limit++;
        
        $new_cum_area = int($cumulative_area * $global_prob_resolution);
        for(my $counter = $new_cum_area; $counter > $prev_cum_area; $counter--)
        {
            $ins_values{($counter / $global_prob_resolution)} = $key;
        }
        $prev_cum_area = $new_cum_area;
        $key++;
        $x = $lower_insert_limit; 
    }
    # do the final half notch 
    for($x = $lower_insert_limit; $x < $upper_insert_limit; $x = $x + $step)
    {
        my $height = $multiplier * exp(-1 * ($x - $global_insert) * ($x - $global_insert) / $two_sig_squared);
        $cumulative_area = $cumulative_area + $height * $step; 
    }
    $new_cum_area = int($cumulative_area * $global_prob_resolution);
    for(my $counter = $new_cum_area; $counter >= $prev_cum_area; $counter--)
    {
         $ins_values{($counter / $global_prob_resolution)} = $upper_insert_limit;
    }
    $global_ins_upper_rand = $cumulative_area;

    return %ins_values;
}

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", 
                             "reference|r:s",
                             "num_reads|n:i",
                             "length|l:i",
                             "insert|i:i",
                             "deviation|d:i",
                             "orientation|o:s",
                             "cutoff|c:i"
                           );
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    #if(!exists $options{''} ) { printParamError (""); }

    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #  
    my ($error) = @_;  
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          }, 
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};
    
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to 
    # checkAndRunCommand
    #
    my $ref = shift;
    
    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}


######################################################################
# MISC

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    sammy.pl

=head1 COPYRIGHT

   copyright (C) Michael Imelfort

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

   Stupidly simple paired read simulation

=head1 SYNOPSIS

    sammy.pl -reference <fasta> -num_reads <number_reads>

      -reference -r FASTA        File containing the reference sequence
      -num_reads -n INT          Number of reads to make
      [-length -l INT ]          Read length [default: 100]
      [-insert -i INT ]          The mean insert [default: 500]
      [-deviation -d INT ]       The standard deviation about the mean [default: 50]
      [-orientation -o STRING ]  The orientation of the reads ('same', 'out', 'in') [default: in]
      [-cutoff -c INT ]          Ignore references shorter than this length [default: 1000]
      [-help -h]                 Displays basic usage information
         
=cut

