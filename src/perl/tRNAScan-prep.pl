#!/usr/bin/perl

=head1 NAME

tRNAscan-prep.pl - prepare tRNAscan raw output for mysql upload

=head1 SYNOPSIS

USAGE: tRNAscan-prep.pl
            --input=/path/to/metagene/raw/output
            --outdir=/output/directory
            --database=/path/to/sqlite3/database.db
          [ --log=/path/to/logfile
            --debug=N
          ]

=head1 OPTIONS

B<--input, -i>
    The full path to metagene raw output

B<--database, -b>
    Path to sqlite3 database

B<--outdir, -o>
    output directory location

B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script is used to prepare tRNAScan results for mysql upload.

=head1  INPUT

The input to this is defined using the --input.  This should point
to the tRNAScan-SE raw output file list.  One file per line.

=head1  CONTACT

    Jaysheel D. Bhavsar
    bjaysheel@gmail.com

=cut

use strict;
use warnings;
use DBI;
use Pod::Usage;
use Data::Dumper;
use Fcntl qw( O_TRUNC O_RDONLY O_RDWR O_CREAT);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

BEGIN {
  use Ergatis::Logger;
}

###############################################################################
my %options = ();
my $results = GetOptions (\%options,
                            'input|i=s',
                            'outdir|o=s',
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

###############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

# check if tRNA output file is not empty
# check if the file is empty.
unless(-s $options{input} > 0){
    print STDERR "This file $options{input} seem to be empty nothing therefore nothing to do.";
    $logger->debug("This file $options{input} seem to be empty nothing therefore nothing to do.");
    exit(0);
}

###############################################################################
#### database connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{database}", "", "", { RaiseError => 1 }) or $logger->logdie($DBI::errstr);

my $sth = $dbh->prepare('SELECT id,name FROM sequence WHERE 1');
$sth->execute;

my %sequence_hash;
while ( my $result = $sth->fetchrow_hashref() ) {
	$sequence_hash{$result->{name}} = $result->{id};
}

## tRNAScan-SE score threshold
my $threshold = 20.00;

my $filename = $options{outdir}."/tRNA.txt";

## open handler to read input file.
open (DAT, "<", $options{input}) || die $logger->logdie("Could not open file $options{input}\n");
open (OUT, ">>", $filename) || die $logger->logdie("Could not open file $filename\n");

#loop through input and upload them to db
while (<DAT>){
    unless (/^#/){
        chomp $_;
        my @info = split(/\t/,$_);

        ## check for duplicate entry and threshold cut off.
        if ($info[8] >= $threshold) {
            $info[4] =~ s/\$([a-z])/\u$1/ig;
            $info[4] =~ s/pse/Pseudo/i;
            $info[4] =~ s/und/Undef/i;

            print OUT join("\t",$sequence_hash{trim($info[0])}, trim($info[1]), trim($info[2]), trim($info[3]),
            	   trim($info[4]), trim($info[5]), trim($info[6]), trim($info[7]), trim($info[8]))."\n";
        }
    }
}

close DAT;
close OUT;

exit(0);

###############################################################################
sub check_parameters {
    my @required = qw(input outdir database);

    foreach my $key (@required) {
        unless ($options{$key}) {
            pod2usage({-exitval => 2,  -message => "ERROR: Input $key not defined", -verbose => 1, -output => \*STDERR});
            $logger->logdie("Input option $key not defined, plesae read perldoc $0\n\n");
        }
    }
}

###############################################################################
sub trim {
    my $str = shift;

    $str =~ s/^\s+//;
    $str =~ s/\s+$//;

    return $str;
}
