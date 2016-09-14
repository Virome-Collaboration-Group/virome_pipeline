#!/usr/bin/perl

=head1 NAME

blast-result-prep.pl - prepare blast btabl file for mysql upload

=head1 SYNOPSIS

USAGE: blast-result-prep.pl
            --input=/path/to/blast/btab/output
            --liblist=/library/list/file/from/db-load-library
            --lookupDir=/dir/where/mldbm/lookup/files/are
          [ --log=/path/to/logfile
            --debug=N
          ]

=head1 OPTIONS

B<--input, -i>
    The full path to blast btab output

B<--liblist, -ll>
    Library list file, and output of db-load-library.

B<--lookupDir, -ld>
    Dir where all lookup files are stored.

B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script is used to load blast results in to MySQL database.

=head1  INPUT

The input to this is defined using the --input.  This should point
to the blast btab output

Fields that are load are

  1   query_name
  2   query_length
  3   algorithm
  4   database_name
  5   hit_name
  6   qry_start
  7   qry_end
  8   hit_start
  9  hit_end
  10  percent_identity
  11  percent_similarity
  12  raw_score
  13  bit_score
  14  hit_description
  15  blast_frame
  16  qry_strand (Plus | Minus)
  17  hit_length
  18  e_value

  if UNIREF100P blast result then following are also added

  19  domain
  20  kingdom
  21  phylum
  22  class
  23  order
  24  family
  25  genus
  26  species
  27  organism
  28  fxn_topHit

=head1  CONTACT

    Jaysheel D. Bhavsar
    bjaysheel@gmail.com

=cut

use strict;
use warnings;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use MLDBM 'DB_File';
use Fcntl qw( O_TRUNC O_RDONLY O_RDWR O_CREAT);
use UTILS_V;
use File::Basename;
use Data::Dumper;

BEGIN {
    use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
                          'output_dir|o=s',
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

#### proceed if file size is greater than 0, mostly a check for rRNA blast.
unless(-s $options{input} > 0){
    print STDERR "This file $options{input} seem to be empty nothing therefore nothing to do.";
    $logger->debug("This file $options{input} seem to be empty nothing therefore nothing to do.");
    exit(0);
}

##############################################################################
my $utils = new UTILS_V;
my $max_id = 0;

#### name
my $prev_seq="";
my $curr_seq="";
my $prev_db="";
my $curr_db="";
my @empty_lineage = ("", "", "", "", "", "", "", "", "");

my $filename = fileparse($options{input}, qr/\.[^.]*/);

#### open handler to read input file.
open (DAT, "<", $options{input}) or $logger->logdie("Could not open file $options{input}");
open (OUT, ">", $options{output_dir} . "/" . $filename .".btab") or $logger->logdie("Could not open file ". $options{output_dir} . "/" . $filename .".btab");

#### database connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{database}", "", "", { RaiseError => 1 }) or $logger->logdie($DBI::errstr);

my $sth = $dbh->prepare('SELECT max(id) as max_id FROM blastp');
$sth->execute;

while ( my $result = $sth->fetchrow_hashref() ) {
	if (defined $result->{max_id}) {
		$max_id = $result->{max_id} + 1;
	} else {
		$max_id = 1;
	}
}

my $timestamp = getTimeStamp();

while (<DAT>) {
	next if ($_ =~ /^#/);

	my $topHit = 0;
	my $fxnHit = 0;
	my $line = $_;

    chomp $line;

    my @info = split (/\t/, $line);

    #### check if self blast result
    if ($info[1] ne $info[6]) {
		#### create an in order array to print
		my @to_print_array = $max_id;

		push(@to_print_array, @info[0..19]);

		#### if this is an UNIREF HIT then append functional taxonomy
		if ($#info > 19) {
			push (@to_print_array, @info[20..28]);

			$fxnHit = $info[29];
		} else {
			push (@to_print_array, @empty_lineage);
		}

		#### check if this is the first hit, and set tophit flag
		$curr_seq = $info[1];
		$curr_db = $info[4];

		if (($curr_seq ne $prev_seq) || ($curr_db ne $prev_db)) {
			$topHit = 1;

			$prev_seq = $curr_seq;
			$prev_db = $curr_db;
		}

		print OUT join("\t", @to_print_array, $topHit, $fxnHit, $timestamp, "0000-00-00 00:00:00", 0). "\n";

		$max_id++;
    } #### end check for self blast.
    else {
		$logger->info("WARNING: SELF BLAST :$info[1]---$info[6]");
    }

} #end while loop

close DAT;
close OUT;

$logger->info("BLAST RESULT PREP for $options{sample} complete");

exit(0);

##############################################################################
sub check_parameters {
    ## at least one input type is required
    unless ($options{input} && $options{output_dir} && $options{database}){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
        $logger->logdie("No input defined, plesae read perldoc $0\n\n");
    }
}

###############################################################################
sub getTimeStamp{
	my $ts = `date +"%Y-%m-%d %H:%M:%S"`;

	chomp $ts;

	return $ts;
}
