#!/usr/bin/perl

=head1 NAME
   nt_fasta_check.pl

=head1 SYNOPSIS

    USAGE: nt_fasta_check.pl

=head1 OPTIONS

B<--fasta,-f>
    input fasta file

B<--outdir,-o>
    output directory to store results.

B<--help,-h>
   This help message

=head1  DESCRIPTION
	Takes as input a fasta file, and library prefix/library file
	Prefix each sequence id with PREFIX from library file, and
	check the quality of read bases.

=head1  INPUT
	Fasta file in .fsa, .fa, .fasta or .txt file format
	OUTPUT dir where updated input file is stored (output dir cannot be same as where input file is)

=head1  OUTPUT
	An updated fasta file with each sequenceId prefixed by PREFIX
	A ref file with original and new sequenceId

=head1  CONTACT
	bjaysheel@gmail.com

==head1 EXAMPLE
   nt_fasta_check.pl -i=input.fsa -o=/output_dir
   or
   nt_fasta_check.pl -i=input.fsa -o=/output_dir

=cut


use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Bio::SeqIO;
use LIBInfo;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

BEGIN {
    use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
                          'output_dir|o=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();
#############################################################################
#### DEFINE GLOBAL VAIRABLES.
##############################################################################
#### make sure everything passed was peachy
&check_parameters(\%options);

my $count=0;
my $filebase = fileparse($options{input},  qr/\.[^.]*/);
my ($warn1,$warn2,$warn3,$warn4,$warn5) = (0,0,0,0,0);

my $final_output=$options{output_dir}."/".$filebase.".edited.fsa";
my $ref_file=$options{output_dir}."/".$filebase.".ref";

open(FOUT, ">", $final_output) or $logger->logdie("Cannot open output file $final_output\n");
open(REF, ">", $ref_file) or $logger->logdie("Cannot open ref output file $ref_file\n");

my $inseq = Bio::SeqIO->new(-file   => $options{input},
                            -format => 'fasta' );

while (my $s = $inseq->next_seq) {
	#my $new_name = name_modifier($s->id, $libObject->{prefix});
    my $new_name = name_modifier($s->id);
	my $qc_flag = freq_cal($s->seq, $s->id);
	my $sequence_string = $s->seq;
	my $ATGC = $sequence_string =~ tr/ATGCatgc/ATGCatgc/;
	my $sequence_size = length($sequence_string);
	my $atgc = $ATGC / $sequence_size;

	if ($qc_flag == 1) {
	    print FOUT ">".$new_name."\n".$s->seq."\n";
	}

	$count++;
}


#if(($warn1/$count) > 0.05) {
#	print STDERR "ERROR 255: Number of ATCG minor warnings exceeds 5\% (".($warn1/$count).")\n";
#	exit(255);
#}

#if(($warn2/$count) > 0.03) {
#	print STDERR "ERROR 256: Number of ATCG major warnings exceeds 3\% (".($warn2/$count).")\n";
#    exit(256);
#}
#if(($warn3/$count) > 0.05) {
#	print STDERR "ERROR 257: Number of N minor warnings exceeds 5\% (".($warn3/$count).")\n";
#    exit(257);
#}
#if(($warn4/$count) > 0.03) {
#	print STDERR "ERROR 258: Number of N major warnings exceeds 3\% (".($warn4/$count).")\n";
#    exit(258);
#}

close(FOUT);
close(REF);

exit(0);
###############################################################################
####  SUBS
###############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(input output_dir);

	foreach my $key (@required) {
		unless ($options{$key}) {
			pod2usage({-exitval => 2,  -message => "ERROR: Input $key not defined", -verbose => 1, -output => \*STDERR});
			$logger->logdie("Inputs not defined, plesae read perldoc $0\n");
		}
	}
}

###############################################################################
sub name_modifier{
	my $id = shift;
	#my $prefix = shift;

    my $name = $id;

	## replace any underscore with dash to identify
	## prefix from rest of the sequence.
        #$id =~ s/_/-/g;

    #### get rid of pipes and slashes in the header too
	$id =~ s/\|/-/g;
	$id =~ s/\//-/g;

    #### get rid of parentheses as well
	$id =~ s/\(//g;
	$id =~ s/\)//g;

    #my $name = $prefix . '_'. $id;

	#### add old id to new id mapping
	print REF $id."\t".$name."\n";
	if (length($id) > 255) {
	    $logger->logdie("ERROR: There was a FASTA header longer than 255 characterss. Here's the culprit: \n $id");
	}

	return $id;
}

###############################################################################
sub freq_cal
{
    my $seq = shift;
    my $seq_name = shift;
    my $flag = 1;

    my $len = length($seq);
    $seq =~ s/X//gi;
    if ($seq =~ m/[^ATCGNRYSWKMBDHV]/i) {
    	# print STDERR "Invalid base(s) in $seq_name";
    	# exit(259);
    	print STDERR "Warning (Major) for INVALID BASES for seq id: $seq_name\n";
    	$flag = 0;
    	$warn5++;
    }

    my $ATCGcount = () = $seq =~ /[ATCG]/ig;
    my $Ncount = () = $seq =~ /[N]/ig;

    my $freq_atcg = ($ATCGcount/$len)*100;
    my $freq_n = ($Ncount/$len)*100;

    if($freq_atcg < 97) {
    	print STDERR "Warning (Minor) for ATCG Frequency (".$freq_atcg."\%) for seq id: $seq_name\n";
    	$warn1++;
    }
    if($freq_atcg < 93) {
    	print STDERR "Warning (Major) for ATCG Frequency (".$freq_atcg."\%) for seq id: $seq_name\n";
    	$flag = 0;
    	$warn2++;
    }
    if($freq_n > 5) {
    	print STDERR "Warning (Major) for N frequency (".$freq_n."\%)  for seq id: $seq_name\n";
    	$flag = 0;
    	$warn4++;
    }
    if($freq_n > 2) {
    	print STDERR "Warning (Minor) for N frequncy (".$freq_n."\%) for seq id: $seq_name\n";
    	$warn3++;
    }
    return $flag;
}
