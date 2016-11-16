#!/usr/bin/perl

=head1 NAME
   libraryHistorgram.pl

=head1 SYNOPSIS

    USAGE: libraryHistogram.pl --server server-name --env dbi [--library libraryId]

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
    Create XML document that contaions information to draw all Fxnal
    count chart. ie: number of hits in keeg, seed, uniref100p, cog, alcame

=head1  INPUT
    The input is defined with --server,  --library.

=head1  OUTPUT
   Updated blastp tables for all/specifed library.

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   libraryHistogram.pl --server calliope --env dbi --library 31

=cut

use strict;
use warnings;
use IO::File;
use DBI;
use LIBInfo;
use UTILS_V;
use XML::Writer;
use Pod::Usage;
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
my $this;
my $cmd = "";

$this->{output_dir} = "$options{outdir}/xDocs/";

## make sure everything passed was peachy
&check_parameters(\%options);

#init db connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{input}", "", "", { RaiseError => 1}) or die $DBI::errstr;

my $libId = 1;

my $blst_sel = $dbh->prepare(q{SELECT DISTINCT b.sequenceId
							 FROM blastp b LEFT JOIN sequence s on s.id=b.sequenceId
					WHERE 	b.e_value < 0.001
						and s.libraryId=?
						and b.database_name = ?
						and b.fxn_topHit=1
					and b.deleted=0 and s.deleted=0
			       ORDER BY sequenceId});

my $xml_file = "FXNAL_OVERVIEW_XMLDOC_".$libId.".xml";
my $id_file = "FXNAL_OVERVIEW_IDDOC_".$libId.".xml";
my $tag_num = 1;

my $xml_out = new IO::File(">". $this->{output_dir} ."/". $xml_file) or $logger->logdie("Could not open file ".$this->{output_dir} ."/". $xml_file." to write");
my $id_out = new IO::File(">". $this->{output_dir} ."/". $id_file) or $logger->logdie("Could not open file ".$this->{output_dir} ."/".$id_file." to write");

my $xml_writer = new XML::Writer(OUTPUT=>$xml_out);
my $id_writer = new XML::Writer(OUTPUT=>$id_out);

$xml_writer->xmlDecl("UTF-8");
$id_writer->xmlDecl("UTF-8");

$xml_writer->startTag("root");
$id_writer->startTag("root");

($xml_writer, $id_writer, $tag_num) = getNodeFor('ACLAME', $libId, $xml_writer, $id_writer, $id_file, $tag_num);
($xml_writer, $id_writer, $tag_num) = getNodeFor('COG', $libId, $xml_writer, $id_writer, $id_file, $tag_num);
($xml_writer, $id_writer, $tag_num) = getNodeFor('UNIREF100P', $libId, $xml_writer, $id_writer, $id_file, $tag_num);
($xml_writer, $id_writer, $tag_num) = getNodeFor('KEGG', $libId, $xml_writer, $id_writer, $id_file, $tag_num);
($xml_writer, $id_writer, $tag_num) = getNodeFor('SEED', $libId, $xml_writer, $id_writer, $id_file, $tag_num);
($xml_writer, $id_writer, $tag_num) = getNodeFor('PHGSEED', $libId, $xml_writer, $id_writer, $id_file, $tag_num);

$xml_writer->endTag("root");
$id_writer->endTag("root");

$xml_writer->end();
$id_writer->end();

$xml_out->close();
$id_out->close();

$logger->info("Functional binning all database for $options{sample} complete");
exit(0);

###############################################################################
####  SUBS
###############################################################################
sub check_parameters {
    my $options = shift;

    my @required = qw(input outdir);

    foreach my $key (@required) {
        unless ($options{$key}) {
            pod2usage({-exitval => 2,  -message => "ERROR: Input $key not defined", -verbose => 1, -output => \*STDERR});
            $logger->logdie("No input defined, plesae read perldoc $0\n\n");
        }
    }

    system ("mkdir -p $options{outdir}/idFiles");
    system ("mkdir -p $options{outdir}/xDocs");

}

###############################################################################
sub timer {
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
    print "Time now: " . $theTime."\n";
}

###############################################################################
sub getNodeFor {
    my ($db, $lib, $xw, $iw, $fname, $tag) = @_;

    $blst_sel->execute($lib,$db);
    my $rslt = $blst_sel->fetchall_arrayref({});

    if (@{$rslt} > 0){
	my $idList = '';
	my $count = 0+@{$rslt};

	foreach my $row (@$rslt) {
	    $idList .= $row->{sequenceId} . ",";
	}

	$idList =~ s/,$//;

	if ($db eq 'UNIREF100P'){
	    $db = 'GO';
	}

	$xw->emptyTag("CATEGORY", 'LABEL'=>$db, 'VALUE'=>$count, 'TAG'=>'TAG_'.$tag, 'IDFNAME'=>$fname);
	$iw->emptyTag("TAG_".$tag, 'IDLIST'=>$idList);
    }

    return $xw, $iw, $tag+1;
}
