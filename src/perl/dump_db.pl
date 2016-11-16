=head1 NAME

dump_db.pl - Dumps the contents of a database into a directory full of table tab files

=head1 SYNOPSIS

USAGE: dump_db.pl
            --info=/Path/to/library.txt
            --outdir=/Path/to/outdir
            --mgol="MGOL_VERSION"
            --uniref="UNIREF_VERSION"
            --pipeline="PIPELINE_VERSION"

=head1 OPTIONS

B<--info,-i>
    The library info file

B<--outdir,-o>
    The output directory

B<--mgol,-m>
    The MgOl BLAST DB version

B<--uniref,-u>
    The UniRef BLAST DB version

B<--pipeline,-p>
    The pipeline version

B<--help,-h>
    This help message

=head1  DESCRIPTION

Dumps the contents of a database into a directory full of table tab files

=head1  INPUT
Output of db-load-library. Essentially a tab-dleimmited file containing:
library_id    library_name    prefix    server    processing_server

=head1  OUTPUT

Directory filled with tab-delimmited table dumps.

=head1  CONTACT

    Daniel Nasko
    dan.nasko@gmail.com

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
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
						  'mgol|m=s',
                          'uniref|u=s',
                          'pipeline|p=s',
                          'pipelineid|r=s',
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
my $outdir = "/opt/output". $filename ."_". timestamp();

system ("mkdir -p $outdir")
system ("mkdir -p $outdir/idFiles");
system ("mkdir -p $outdir/xDocs");
###############################################################################
$logger->info("Database dump for $filename started");

my $dbh = DBI->connect("dbi:SQLite:dbname=$options{database}", "", "", { RaiseError => 1}) or die $DBI::errstr;

###################################
## Let's Gather Some Information ##
###################################
my $library_id = 1;
my $library_name = $filename;
my $prefix = "NAN";

my @tables = ('blastp','blastn','sequence','statistics','sequence_relationship','tRNA');

foreach my $table (@tables) {
    print "Dumping table $table\n";

    open(OUT, ">", "$outdir/${table}.dump.tbl") or $logger->logdie("Could not open file to write $outdir/${table}.dump.tbl");

    print OUT "PRAGMA synchronous=OFF;\n";
    print OUT "PRAGMA count_changes=OFF;\n";
    print OUT "PRAGMA journal_mode=MEMORY;\n";
    print OUT "PRAGMA temp_store=MEMORY;\n";
    print OUT ".separator \"\\t\"\n";
    print OUT ".headers on\n";
    print OUT ".out $outdir/$table.tab\n";
    print OUT "select * from $table\n";

    close(OUT);

    $cmd = "sqlite3 $options{database} < $outdir/${table}.dump.tbl";

    system($cmd);

    if (( $? >> 8 ) != 0 ){
    	print STDERR "command failed: $!\n";
    	print STDERR $cmd."\n";
    	exit($?>>8);
    }

    #### add # to header of each table output
    $cmd = "sed -i '1s/^/#/' $outdir/$table.tab";
    system($cmd);
}

$cmd = "cp $options{output_dir}/xDocs/* $outdir/xDocs";
system($cmd);

$cmd = "cp $options{output_dir}/idFiles/* $outdir/idFiles";
system($cmd);

#####################################################################
## Print out the version control info to the version_info.txt file ##
#####################################################################

open(OUT, ">", "$outdir/version_info.txt") || die "\n Cannot open the file: $outdir/version_info.txt\n";
print OUT "fxndbLookupVersion=" . $uniref . "\n";
print OUT "mgolVersion=" . $mgol . "\n";
print OUT "pipelineVersion=" . $pipeline . "\n";
print OUT "prefix=" . $prefix . "\n";
print OUT "id=" . $library_id . "\n";
close(OUT);

#########################
## Create the Tar Ball ##
#########################
if (-e "$outdir.tar.gz" ) {
    print `rm $outdir.tar.gz`;
}
print `tar -czvf $outdir.tar.gz --directory=$outdir $outdir`;

$logger->info("Database dump for $filename completed");
exit(0);

###############################################################################
sub check_parameters {
    my @required = qw(input mgol uniref pipeline pipelineid database output_dir);

    foreach my $key (@required) {
        unless ($options{$key}) {
            pod2usage({-exitval => 2,  -message => "ERROR: Input $key not defined", -verbose => 1, -output => \*STDERR});
            $logger->logdie("No input defined, plesae read perldoc $0\n\n");
        }
    }
}

###############################################################################
sub timestamp {
	my @months   = qw(01 02 03 04 05 06 07 08 09 10 11 12);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my (
		$second,     $minute,    $hour,
		$dayOfMonth, $month,     $yearOffset,
		$dayOfWeek,  $dayOfYear, $daylightSavings
	  )
	  = localtime();

	my $year    = 1900 + $yearOffset;
	my $theTime = $year ."". $months{$month} ."". $dayOfMonth ."". $hour ."". $minute ."". $second;

    return $theTime;
}
