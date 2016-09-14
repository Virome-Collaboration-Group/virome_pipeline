#!/usr/bin/perl

=head1 NAME

db-upload.pl: Uplaod file into mysql db

=head1 SYNOPSIS

USAGE: db-upload.pl
            --input=/file/preped/for/mysqlimport
			--table=/tablename
			--env=/env/where/executing
			--outdir=/output/dir/loc
			[ --log=/path/to/logfile
			--debug=N]

=head1 OPTIONS

B<--input, -i>
    The full path to tab delimited file prepared for import.

B<--table, -t>
    mysql db table name

B<--database, -b>
    location of sqlite3 database file

B<--outdir, -o>
    output directory to write sqlite import file

B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script is used to upload info into db.

=head1  INPUT

=head1  CONTACT

    Jaysheel D. Bhavsar
    bjaysheel@gmail.com

=cut

use strict;
use warnings;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use UTILS_V;

BEGIN {
  use Ergatis::Logger;
}

##############################################################################
my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output_dir|o=s',
						  'table|t=s',
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

my $filename = fileparse($options{input}, qr/\.[^.]*/);
###############################################################################
$logger->info("Database upload for $filename in $options{table} started");

my $cmd = "";
my $this;

open(OUT, ">", "$options{output_dir}/$options{table}") or $logger->logdie("Could not open file to write $options{output_dir}/$options{table}");

print OUT "PRAGMA synchronous=OFF;\n";
print OUT "PRAGMA count_changes=OFF;\n";
print OUT "PRAGMA journal_mode=MEMORY;\n";
print OUT "PRAGMA temp_store=MEMORY;\n";
print OUT ".separator \"\\t\"\n";
print OUT ".import $options{input} $options{table}\n";

close(OUT);

$cmd = "sqlite3 $options{database} < $options{output_dir}/$options{table}";

execute_cmd($cmd);

if (( $? >> 8 ) != 0 ){
	print STDERR "command failed: $!\n";
	print STDERR $cmd."\n";
	exit($?>>8);
}

$logger->info("Database upload for $filename in $options{table} completed");
exit(0);

###############################################################################
sub check_parameters {
  ## at least one input type is required
	unless ( $options{input} && $options{table} && $options{database} && $options{outdir}) {
		pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
		$logger->logdie("No input defined, plesae read perldoc $0\n\n");
		exit(1);
	}
}
