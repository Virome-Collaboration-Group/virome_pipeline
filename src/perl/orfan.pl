#!/usr/bin/perl

=head1 NAME
   db-load-nohit.pl

=head1 SYNOPSIS

    USAGE: db-load-nohit.pl --server server-name --env dbi [--library libraryId]

=head1 OPTIONS

B<--server,-s>
   Server name from where MGOL blastp records are updated

B<--library,-l>
    Specific libraryId whoes MGOL hit_names to updates

B<--env,-e>
    Specific environment where this script is executed.  Based on these values
    db connection and file locations are set.  Possible values are
    igs, dbi, ageek or test

B<--help,-h>
   This help message

=head1  DESCRIPTION
    Add no-hit values to DB

=head1  INPUT
    The input is defined with --server, --library.

=head1  OUTPUT
   Add no-hit infor to blast tables.

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   no-hit.pl --server calliope --env dbi --library 31

=cut

se strict;
use warnings;
use DBI;
use Switch;
use LIBInfo;
use Pod::Usage;
use Data::Dumper;
use UTILS_V;
use MLDBM 'DB_File';
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

BEGIN {
  use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
            			  'input|i=s',
            			  'outdir|o=s',
            			  'log|l=s',
            			  'debug|d=s',
            			  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ) {
  pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				 'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

$logger->info("ORFAN library stats for $options{sample} started");

##############################################################################
#### DEFINE GLOBAL VAIRABLES.
##############################################################################
my $this;
my $cmd = "";
my $max_id = -1;

$this->{output_dir} = "$options{outdir}/orfan/";
$this->{input_dir} = "$options{outdir}/idFiles/";

## make sure everything passed was peachy
&check_parameters(\%options);

#init db connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{input}", "", "", { RaiseError => 1}) or die $DBI::errstr;


my $filename = "$this->{output_dir}/all.tab";
my $libId = 1;
my $seq_info = $dbh->prepare(qq{SELECT name,size,id FROM sequence where id=?});

my $sth = $dbh->prepare('SELECT max(id) as max_id FROM blastp');
$sth->execute;

while ( my $result = $sth->fetchrow_hashref() ) {
	if (defined $result->{max_id}) {
		$max_id = $result->{max_id} + 1;
	} else {
		$max_id = 1;
	}
}

open( OUT, ">", $filename ) or die $logger->logdie("Could not open file $filename");

####open orf file and get all ids.
open( ORF, "<", $this->{input_dir}."/orfanList_".$libId.".txt" ) or $logger->logdie("Could not open file " .$this->{output_dir}. "/idFiles/orfanList_" .$libId. ".txt");
my $orfs = <ORF>;
my @sequenceIDs = split( /,/, $orfs );

####loop through all orfan ids and get seqeunce info to insert into blast table.
foreach my $id (@sequenceIDs) {
	$seq_info->execute($id);
	my ( $name, $size, $sequenceId ) = $seq_info->fetchrow_array();

	####output data to tab file for import.
	print OUT join("\t", $max_id, $sequenceId,
					   $name, $size, "BLASTP", "NOHIT", 0,
					   "", "Sequence has no homologs in KNOWN or ENVIRONMENTAL database",
					   "", "", "", "", 0.0, 0.0,
					   0, 0, "", "", "", 0.001,
					   "", "", "", "", "", "", "", "", "",
					   1, 0, getTimeStamp(), "0000-00-00 00:00:00", 0, 0);
	print OUT "\n";

	$max_id++;
}

close(OUT);

#### move updated sqlite database back from local to NFS mount
#$cmd = "mv $this->{scratch}/$config->{ToolInfo}->{database_name}->{value} $this->{output_dir}/$config->{ToolInfo}->{database_name}->{value}";
#execute_cmd($cmd);


$logger->info("ORFAN library stats for $options{sample} completed");
exit(0);

############################################################################
sub check_parameters {
    my $options = shift;

    my @required = qw(input output);

    foreach my $key (@required) {
        unless ($options{$key}) {
            pod2usage({-exitval => 2,  -message => "ERROR: Input $key not defined", -verbose => 1, -output => \*STDERR});
            $logger->logdie("No input defined, plesae read perldoc $0\n\n");
        }
    }

    system ("mkdir -p $options{outdir}/idFiles");
    system ("mkdir -p $options{outdir}/xDocs");

}

#############################################################################
sub create_dir_struct {
	my $options = shift;

	my $dir = $this->{output_dir};
	my $str = "";
	my @arr = split(/\//, $dir);

	#### assuming path always starts with /
	#### in which case first array element is empty
	@arr = @arr[1..$#arr];

	foreach my $d (@arr){
		$str .= "/$d";

		unless (-d $str) {
			mkdir $str, 0770 or $logger->logdie("Could not create dir $str");
			$logger->info("mkdir $str, 0770");
		}
	}

	$dir = $this->{scratch};
	$str = "";
	@arr = split(/\//, $dir);

	#### assuming path always starts with /
	#### in which case first array element is empty
	@arr = @arr[1..$#arr];

	foreach my $d (@arr){
		$str .= "/$d";

		unless (-d $str) {
			mkdir $str, 0770 or $logger->logdie("Could not create dir $str");
			$logger->info("mkdir $str, 0770");
		}
	}
}

#############################################################################
sub execute_cmd {
	my $cmd = shift;

	$logger->info("$cmd");
	system($cmd);

	my $c=0;
	while (( $? >> 8 ) != 0 ){

		#### work around for file system sync lag
		#### busy wait and try the command few times before failing.
		if ($c < 5) {
			$c++;
			sleep 5;
			system($cmd);
		} else {
			$logger->logdie("ERROR: Following command failed to execute. Exiting execution of workflow\n$cmd");
			exit(100);
		}
	}
}

###############################################################################
sub getTimeStamp{
	my $ts = `date +"%Y-%m-%d %H:%M:%S"`;

	chomp $ts;

	return $ts;
}
