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
use File::Path qw(make_path remove_tree mkpath);

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
                          'pstore|s=s',
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

my $persistent_outdir = $options{pstore};
my $dirname = fileparse($persistent_outdir, qr/\.[^.]*/);

my $rindex = rindex($dirname, "_");
my $filename = substr($dirname, 0, $rindex);  #### remove timestamp

my $cmd = "";

#### create idFiles and xDocs dir
make_path($persistent_outdir."/idFiles");
make_path($persistent_outdir."/xDocs");
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

    open(OUT, ">", "$persistent_outdir/${table}.dump.tbl") or $logger->logdie("Could not open file to write $persistent_outdir/${table}.dump.tbl");

    print OUT "PRAGMA synchronous=OFF;\n";
    print OUT "PRAGMA count_changes=OFF;\n";
    print OUT "PRAGMA journal_mode=MEMORY;\n";
    print OUT "PRAGMA temp_store=MEMORY;\n";
    print OUT ".separator \"\\t\"\n";
    print OUT ".headers on\n";
    print OUT ".out $persistent_outdir/$table.tab\n";
    print OUT "select * from $table;\n";

    close(OUT);

    $cmd = "sqlite3 $options{database} < $persistent_outdir/${table}.dump.tbl";

    system($cmd);

    if (( $? >> 8 ) != 0 ){
    	print STDERR "command failed: $!\n";
    	print STDERR $cmd."\n";
    	exit($?>>8);
    }

    #### add # to header of each table output
    $cmd = "sed -i '1s/^/#/' $persistent_outdir/$table.tab";
    system($cmd);
}

$cmd = "cp $options{output_dir}/xDocs/* $persistent_outdir/xDocs";
system($cmd);

$cmd = "cp $options{output_dir}/idFiles/* $persistent_outdir/idFiles";
system($cmd);

#####################################################################
## Print out the version control info to the version_info.txt file ##
#####################################################################

open(OUT, ">", "$persistent_outdir/version_info.txt") || die "\n Cannot open the file: $persistent_outdir/version_info.txt\n";
print OUT "fxndbLookupVersion=" . $options{uniref} . "\n";
print OUT "mgolVersion=" . $options{mgol} . "\n";
print OUT "pipelineVersion=" . $options{pipeline} . "\n";
print OUT "prefix=" . $prefix . "\n";
print OUT "id=" . $library_id . "\n";
close(OUT);

#### remove all .tbl files
$cmd = "rm -rf $persistent_outdir/*.dump.tbl";
system($cmd);

#########################
## Create the Tar Ball ##
#########################
#### persistent_outdir contatins full path plus dirname,
#### final output is a tar file stored in persistent_outdir plus exe dir name
#### and extention of .tar
if (-e "$persistent_outdir.tar" ) {
    $cmd = "rm $persistent_outdir.tar";
    system($cmd);
}

#### create tar without /opt/output in the path.
$cmd = "tar --exclude=\"$dirname/processing.sqlite3\"";
$cmd .= " -czvf $persistent_outdir.tar.gz -C /opt/output $dirname";
system($cmd);

#### get md5sum and touch a file with that name.
my $md5sum = `md5sum $persistent_outdir.tar.gz`;
chomp $md5sum;

#### extract only the checksum
$md5sum =~ s/ $persistent_outdir.tar.gz//;

#### trim any spaces there might be;
$md5sum =~ s/^\s+//;
$md5sum =~ s/\s+$//;

#### touch a file with md5sum.
$cmd = "touch /opt/output/$md5sum";
system($cmd);

#### create a tarball for md5sum file and tar.gz file.
$cmd = "tar -cvf $persistent_outdir.tar -C /opt/output $dirname.tar.gz $md5sum";
system($cmd);

#### remove unwanted files;
unless ($options{debug}){
    $cmd = "rm -rf $persistent_outdir.tar.gz /opt/output/$md5sum $persistent_outdir";
    system($cmd);
}

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

    $options{debug} = 0 unless (defined $options{debug});
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
	my $theTime = $year ."". $months[$month] ."". $dayOfMonth ."". $hour ."". $minute ."". $second;

    return $theTime;
}
