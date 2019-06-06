#!/usr/bin/perl

=head1 NAME
   viromeClassification.pl

=head1 SYNOPSIS

    USAGE: viromeClassification.pl -server server-name --env dbi [--library libraryId]

=head1 OPTIONS

B<--server,-s>
   Server name from where MGOL blastp records are updated

B<--library,-l>
    Specific libraryId whoes MGOL hit_names to updates

B<--help,-h>
   This help message

=head1  DESCRIPTION
    Create XML document that contaions information to draw virome classification
    and database breakdown charts

=head1  INPUT
    The input is defined with --server,  --library.

=head1  OUTPUT
   Updated blastp tables for all/specifed library.

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   viromeClassification.pl --server calliope --env dbi --library 31

=cut

use strict;
use warnings;
use IO::File;
use POSIX qw/ceil/;
use DBI;
use UTILS_V;
use XML::Writer;
use Pod::Usage;
use File::Path qw(make_path remove_tree mkpath);
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
#############################################################################
my $this;

$this->{output_dir} = "$options{outdir}/xDocs/";
$this->{input_dir} = "$options{outdir}/idFiles/";

## make sure everything passed was peachy
&check_parameters(\%options);

#init db connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{input}", "", "", { RaiseError => 1}) or die $DBI::errstr;

my $libId = 1;

##############################################################################

my $stat_sel = $dbh->prepare(q{SELECT	s.read_cnt,s.read_mb,s.complete_cnt,
					s.complete_mb,s.complete_id,s.incomplete_cnt,
					s.incomplete_mb,s.incomplete_id,s.lackstop_cnt,
					s.lackstop_mb,s.lackstop_id,s.lackstart_cnt,
					s.lackstart_mb,s.lackstart_id,s.archaea_cnt,
					s.archaea_mb,s.archaea_id,s.bacteria_cnt,
					s.bacteria_mb,s.bacteria_id,s.phage_cnt,
					s.phage_mb,s.phage_id,s.tRNA_cnt,
					s.tRNA_id,s.rRNA_cnt,s.rRNA_id,s.orfan_cnt,
					s.orfan_id,s.allviral_cnt,s.allviral_id,
					s.topviral_cnt,s.topviral_id,
					s.allmicrobial_cnt,s.allmicrobial_id,
					s.topmicrobial_cnt,s.topmicrobial_id,
					s.fxn_cnt,s.fxn_id,
					s.unassignfxn_cnt,s.unassignfxn_id,
					s.libraryId
				FROM	statistics s
				WHERE	s.libraryId = ?});

my $mgol_sel = $dbh->prepare(q{SELECT	distinct (b.sequenceId), b.hit_name
				FROM	blastp b
					INNER JOIN
					sequence s on s.id=b.sequenceId
				WHERE	b.e_value < 0.001
					and s.libraryId = ?
					and b.sys_topHit=1
					and b.database_name = 'METAGENOMES'
				ORDER BY b.sequenceId, b.hit_name});

my $xml_file = "VIRClass_XMLDOC_".$libId.".xml";
my $db_file = "DBBreakdown_XMLDOC_".$libId.".xml";
my $db_id_file = "DBBreakdown_IDDOC_".$libId.".xml";

my $xml_out = new IO::File(">".$this->{output_dir}."/".$xml_file) or $logger->logdie("Could not open file ".$this->{output_dir}."/".$xml_file." to write");

my $xml_writer = new XML::Writer(OUTPUT=>$xml_out);
$xml_writer->xmlDecl("UTF-8");
$xml_writer->startTag("root");

$stat_sel->execute($libId);
my $rslt = $stat_sel->fetchall_arrayref({});
my $fxnStruct;

foreach my $rec (@$rslt){
	my $tag=1;
	my $total = $rec->{'tRNA_cnt'} + $rec->{'rRNA_cnt'} + $rec->{'fxn_cnt'} +
		$rec->{'unassignfxn_cnt'} + $rec->{'topviral_cnt'} +
		$rec->{'allviral_cnt'} + $rec->{'topmicrobial_cnt'} +
		$rec->{'allmicrobial_cnt'} + $rec->{'orfan_cnt'};

	$xml_writer->emptyTag("CATEGORY", 'LABEL'=>'tRNA',
				'CAT'=>'tRNA',
				'VALUE'=>($rec->{'tRNA_cnt'} > 0) ? ceil(($rec->{'tRNA_cnt'}/$total)*100) : 0);

	$xml_writer->emptyTag("CATEGORY", 'LABEL'=>'rRNA',
				'CAT'=>'rRNA',
				'VALUE'=>($rec->{'rRNA_cnt'} > 0) ? ceil(($rec->{'rRNA_cnt'}/$total)*100) : 0);

	$xml_writer->emptyTag("CATEGORY", 'LABEL'=>'Possible functional protein',
				'CAT'=>'fxn',
				'VALUE'=>($rec->{'fxn_cnt'} > 0) ? ceil(($rec->{'fxn_cnt'}/$total)*100) : 0);

	$xml_writer->emptyTag("CATEGORY", 'LABEL'=>'Unassignfxn protein',
				'CAT'=>'unassignfxn',
				'VALUE'=>($rec->{'unassignfxn_cnt'} > 0) ? ceil(($rec->{'unassignfxn_cnt'}/$total)*100) : 0);

	$xml_writer->emptyTag("CATEGORY", 'LABEL'=>'Top-hit viral',
				'CAT'=>'topviral',
				'VALUE'=>($rec->{'topviral_cnt'} > 0) ? ceil(($rec->{'topviral_cnt'}/$total)*100) : 0);

	$xml_writer->emptyTag("CATEGORY", 'LABEL'=>'Viral only',
				'CAT'=>'allviral',
				'VALUE'=>($rec->{'allviral_cnt'} > 0) ? ceil(($rec->{'allviral_cnt'}/$total)*100) : 0);

	$xml_writer->emptyTag("CATEGORY", 'LABEL'=>'Top-hit microbial',
				'CAT'=>'topmicrobial',
				'VALUE'=>($rec->{'topmicrobial_cnt'} > 0) ? ceil(($rec->{'topmicrobial_cnt'}/$total)*100) : 0);

	$xml_writer->emptyTag("CATEGORY", 'LABEL'=>'Microbial only',
				'CAT'=>'allmicrobial',
				'VALUE'=>($rec->{'allmicrobial_cnt'} > 0) ? ceil(($rec->{'allmicrobial_cnt'}/$total)*100) : 0);

	$xml_writer->emptyTag("CATEGORY", 'LABEL'=>'ORFan',
				'CAT'=>'orfan',
				'VALUE'=>($rec->{'orfan_cnt'} > 0) ? ceil(($rec->{'orfan_cnt'}/$total)*100) : 0);
}

my $bothList = '';
my $mgolList = '';
my $uniList = '';
my $mgolCount = 0;
my $bothCount = 0;
my $uniCount = 0;

#### get uniref/fxnal ids.
open (DAT, "<", $this->{input_dir}."/fxnIdList_".$libId.".txt") or $logger->logdie("Could not open file $this->{output_dir}/fxnIdList_${libId}.txt");
my @data = <DAT>;

my @fxnIds = split(/,/, $data[0]);

#### add fxn ids to struct.
foreach my $id_1 (@fxnIds){
	if (!exists $fxnStruct->{$id_1}){
		$fxnStruct->{$id_1} = 1;
		$uniCount ++;
		$uniList .= $id_1.',';
	}
}
close(DAT);

#### get unclassified id's and add
open (DAT, "<", $this->{input_dir}."/unClassIdList_".$libId.".txt") or $logger->logdie("Could not open file $this->{output_dir}/unClassIdList_${libId}.txt");
undef @data;
@data = <DAT>;

my @unassIds = split(/,/, $data[0]);

#### add unassigned ids to struct.
foreach my $id_2 (@unassIds){
	if (!exists $fxnStruct->{$id_2}){
		$fxnStruct->{$id_2} = 1;
		$uniCount ++;
		$uniList .= $id_2.',';
	}
}
close(DAT);

#### get orfan ids
open (DAT, "<", $this->{input_dir}."/orfanList_".$libId.".txt") or $logger->logdie("Could not open file $this->{output_dir}/orfanList_${libId}.txt");
undef @data;
@data = <DAT>;

my @orfIds = split(/,/, $data[0]);
my $orfList = $data[0];
close(DAT);

#### get mgol ids. and sort out which seq exist in both uniref and mgol
$mgol_sel->execute($libId);
$rslt = $mgol_sel->fetchall_arrayref({});

foreach my $rec (@$rslt){
	$mgolList .= $rec->{'sequenceId'}.",";
	$mgolCount ++;
	if (defined $fxnStruct->{$rec->{'sequenceId'}}){
		$bothList .= $rec->{'sequenceId'}.",";
		$bothCount ++;
	}
}

$bothList =~ s/,$//;
$mgolList =~ s/,$//;
$uniList =~ s/,$//;

#### write out xml for db breakdown.
my $db_out = new IO::File(">". $this->{output_dir}."/".$db_file) or $logger->logdie("Could not open file $this->{output_dir}/$db_file to write");

$xml_writer->endTag("root");
$xml_writer->end();
$xml_out->close();

my $db_writer = new XML::Writer(OUTPUT=>$db_out);
$db_writer->xmlDecl("UTF-8");
$db_writer->startTag("root");

$db_writer->emptyTag('CATEGORY', 'LABEL'=>'Uniref100P',
				 'VALUE'=>$uniCount,
				 'TAG'=>'TAG_1',
				 'IDFNAME'=>$db_id_file);
$db_writer->emptyTag('CATEGORY', 'LABEL'=>'Metagenome Online',
				 'VALUE'=>$mgolCount,
				 'TAG'=>'TAG_2',
				 'IDFNAME'=>$db_id_file);
$db_writer->emptyTag('CATEGORY', 'LABEL'=>'Both',
				 'VALUE'=>$bothCount,
				 'TAG'=>'TAG_3',
				 'IDFNAME'=>$db_id_file);
$db_writer->emptyTag('CATEGORY', 'LABEL'=>'Orfans',
				 'VALUE'=>$#orfIds,
				 'TAG'=>'TAG_4',
				 'IDFNAME'=>$db_id_file);

#write out db ids breakdown.
my $db_id_out = new IO::File(">". "$this->{output_dir}/$db_id_file") or $logger->logdie("Could not open file $this->{output_dir}/$db_id_file to write");

$db_writer->endTag("root");
$db_writer->end();
$db_out->close();

my $db_id_writer = new XML::Writer(OUTPUT=>$db_id_out);
$db_id_writer->xmlDecl("UTF-8");
$db_id_writer->startTag("root");

$db_id_writer->emptyTag('TAG_1', 'IDLIST'=>$uniList);
$db_id_writer->emptyTag('TAG_2', 'IDLIST'=>$mgolList);
$db_id_writer->emptyTag('TAG_3', 'IDLIST'=>$bothList);
$db_id_writer->emptyTag('TAG_4', 'IDLIST'=>$orfList);

$db_id_writer->endTag('root');
$db_id_writer->end();
$db_id_out->close();

$logger->info("Sequence classification by VIROME categories started");
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

    make_path($options{outdir}."/idFiles");
    make_path($options{outdir}."/xDocs");
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
