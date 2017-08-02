#!/usr/bin/perl

=head1 NAME

sequence-prep.pl - prepare sequence info for upload to db

=head1 SYNOPSIS

USAGE: sequence-prep.pl
            --input=/path/to/fasta
			--outdir=/output/dir
			--libListFile=/library/list/file
			--type=sequence type
            --outdir=/output/dir
          [ --log=/path/to/logfile
            --debug=N ]

=head1 OPTIONS

B<--input, -i>
    The full path to fasta sequence file.
    # start a comment.
    # File format
    >ABC125234 ....
    # where ABC is three letter library prefix.

B<--outdir, -od>
    Output dir where sequence prep file is uploaded

B<--libListFile, -ll>
    Library list file eg. from db-load-library.

B<--type, -t>
    sequence type read=1, rRNA=2, orf (aa)=3, orf (dna)=4

B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script is used to prepare sequence for mysql upload.

=head1  INPUT

The input to this is defined using the --input.  This should point
to the fasta file containing sequence(s).  Input must be a muilt fasta
file, each file containg sequences for one library.

=head1  CONTACT

    Jaysheel D. Bhavsar
    bjaysheel@gmail.com

=cut

use strict;
use warnings;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Bio::SeqIO;

BEGIN {
  use Ergatis::Logger;
}

##############################################################################
my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'outdir|o=s',
						  'typeId|t=s',
                          'database|b=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
##############################################################################

## make sure everything passed was peachy
&check_parameters(\%options);

unless(-s $options{input} > 0){
	print STDERR "This file $options{input} seem to be empty nothing therefore nothing to do.";
	$logger->debug("This file $options{input} seem to be empty nothing therefore nothing to do.");
	exit(0);
}

##############################################################################

my $filename = $options{outdir}."/sequence.txt";
my $max_id = 0;

#### database connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{database}", "", "", { RaiseError => 1 }) or $logger->logdie($DBI::errstr);

my $sth = $dbh->prepare('SELECT max(id) as max_id FROM sequence');
$sth->execute;

while ( my $result = $sth->fetchrow_hashref() ) {
	if ( (defined $result->{max_id}) && ($result->{max_id} =~ /^[+-]?\d+$/) ) {
		$max_id = $result->{max_id} + 1;
	} else {
		$max_id = 1;
	}
}

open(FHD, "<", $options{input}) || $logger->logdie("Could not open input file $options{input}");
open (OUT, ">>", $filename) || $logger->logdie("Could not open file $filename");

my $header = "";
my $seq = "";

while(<FHD>) {
    chomp $_;

    #### skip coments
    next if ($_ =~ /^#/);

    #### if header
    if ($_ =~ /^>/) {
        #### if we encoutner > (start of new seq) and there is
        #### data in seq then process previous sequence before current one
        #### this situaion will occur after first seq is read in.
        if (length($seq)) {
            my $gc = &calculate_gc($seq);
            my $timestamp = getTimeStamp();
            my $id = $header;
            $id =~ s/\s.*//;
            my $desc = $header;
            $desc =~ s/$id\s//;

        	print OUT join("\t", $max_id, "1", $id, $desc, $gc, $seq, length($seq), "", $timestamp, "0000-00-00 00:00:00", 0, $options{typeId})."\n";

            $max_id++;

            #### reset seq to hold new seq data.
            $seq = "";
        }

        #### read new header
        $header = $_;
        $header =~ s/^>//;
    } else {
        $seq .= $_;
    }
}

#### process last seq that was read before end of while loop
if (length($seq)) {
    my $gc = &calculate_gc($seq);
    my $timestamp = getTimeStamp();
    my $id = $header;
    $id =~ s/\s.*//;
    my $desc = $header;
    $desc =~ s/$id\s//;

    print OUT join("\t", $max_id, "1", $id, $desc, $gc, $seq, length($seq), "", $timestamp, "0000-00-00 00:00:00", 0, $options{typeId})."\n";
}

close(FHD);

#### use Bio::SeqIO to parse and handle fasta file.
#my $fsa = Bio::SeqIO->new( -file   => $options{input},
#						   -format => 'fasta' );

#while (my $seq = $fsa->next_seq) {
#	my $gc = &calculate_gc($seq->seq());
#    my $timestamp = getTimeStamp();

#	print OUT join("\t", $max_id, "1", $seq->id, $seq->desc, $gc, $seq->seq(), $seq->length(), "", $timestamp, "0000-00-00 00:00:00", 0, $options{typeId})."\n";

#    $max_id++;
#}

close OUT;
exit(0);

###############################################################################
sub check_parameters {
    my @required = qw(input outdir typeId);

    foreach my $key (@required) {
        unless ($options{$key}) {
            pod2usage({-exitval => 2,  -message => "ERROR: Input $key not defined", -verbose => 1, -output => \*STDERR});
            $logger->logdie("No input defined, plesae read perldoc $0\n\n");
        }
    }

	if (($options{typeId} <= 0) || ($options{typeId} > 4)){
		pod2usage({-exitval => 2, -message => "Sequence type $options{type} not recognized valid values: 1,2,3,4", -verbose => 1, -output => \*STDERR});
		$logger->logdie("Sequence type $options{type} not recognized valid values: 1,2,3,4\n");
		exit(1);
	}
}

###############################################################################
sub calculate_gc{
	my $bases = shift;

	my $c = $bases;
	my $g = $bases;

	$c =~ s/[ATG]//ig;
	$g =~ s/[ATC]//ig;

	my $gc_percent = ((length($c)+length($g))/length($bases)) * 100;

	return $gc_percent;
}

###############################################################################
sub getTimeStamp{
	my $ts = `date +"%Y-%m-%d %H:%M:%S"`;

	chomp $ts;

	return $ts;
}
