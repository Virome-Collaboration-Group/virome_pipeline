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
    Create XML document that contaions information to draw histogram
    of GC, read sizes and orf sizes to be displayed on VIROME library
    statistics page.

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
## make sure everything passed was peachy
&check_parameters(\%options);

#init db connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{input}", "", "", { RaiseError => 1}) or die $DBI::errstr;
my $this;
my $cmd = "";

$this->{output_dir} = "$options{outdir}/xDocs";

my $libId = 1;

my $orf_sel = $dbh->prepare(q{SELECT s.size*3 as hval
							  FROM 	sequence s
								INNER JOIN
									sequence_relationship sr on s.id = sr.objectId
							  WHERE s.libraryId = ?
								and sr.typeId = 3
								and s.deleted = 0
							  ORDER BY hval desc});

my $read_sel = $dbh->prepare(q{SELECT s.size as hval
							  FROM 	sequence s
								INNER JOIN
									sequence_relationship sr on s.id = sr.objectId
							  WHERE s.libraryId = ?
								and sr.typeId = 1
								and s.deleted = 0
							  ORDER BY hval desc});

my $gc_sel = $dbh->prepare(q{SELECT s.gc as hval
							  FROM 	sequence s
								INNER JOIN
									sequence_relationship sr on s.id = sr.objectId
							  WHERE s.libraryId = ?
								and sr.typeId = 1
								and s.deleted = 0
							  ORDER BY hval desc});
my $rslt;

$orf_sel->execute($libId);
$rslt = $orf_sel->fetchall_arrayref({});
binORFnREADs($rslt, $libId, "ORF");

$read_sel->execute($libId);
$rslt = $read_sel->fetchall_arrayref({});
binORFnREADs($rslt, $libId, "READ");

$gc_sel->execute($libId);
$rslt = $gc_sel->fetchall_arrayref({});
binGC($rslt, $libId, "GC");

$logger->info("Library histogram for $options{sample} completed");

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
sub binORFnREADs {
    my ($rslt, $lib, $type) = @_;

    if (@{$rslt} > 0){
        my $lastIdx = 0 + @{$rslt} -1;

        my $k = 30;
        my $min = $rslt->[$lastIdx]->{hval};
        my $max = $rslt->[0]->{hval}; #not order is desc hense max is first.
        my $bin = int(($max-$min)/$k);

        my @range_arr = ();
        for (my $i=0; $i<$k; $i++){
            push @range_arr, {'bin' => ($min +  ($bin*$i)),
                              'count' => 0};
        }

        my @t_sort = sort { $b->{bin} <=> $a->{bin} } @range_arr;

        foreach my $row (@$rslt) {
            for (my $i=0; $i<=$#t_sort; $i++){
                if ($row->{hval} >= $t_sort[$i]->{'bin'}){
                    $t_sort[$i]->{'count'}++;
                    $i = $#t_sort;
                }
            }
        }

        @range_arr = sort { $a->{bin} <=> $b->{bin} } @t_sort;
        printXML($lib, $type, \@range_arr);

        #foreach $bin (@range_arr){
        #    print $bin->{bin}."\t".$bin->{count}."\n";
        #}
    } else {
        print "Nothing to do, empty result set\n";
    }
}

###############################################################################
sub binGC {
    my ($rslt, $lib, $type) = @_;
    if (@{$rslt} > 0){
        my $lastIdx = 0 + @{$rslt} -1;

        my $k = 21;
        my $min = 0;
        my $max = 100;
        my $bin = 5;

        my @range_arr;
        for (my $i=0; $i<$k; $i++){
            push @range_arr, {'bin' => ($min +  ($bin*$i)),
                              'count' => 0};
        }

        my @t_sort = sort { $b->{bin} <=> $a->{bin} } @range_arr;

        foreach my $row (@$rslt) {
            for (my $i=0; $i<=$#t_sort; $i++){
                if ($row->{hval} >= $t_sort[$i]->{'bin'}){
                    $t_sort[$i]->{'count'}++;
                    $i = $#t_sort;
                }
            }
        }

        @range_arr = sort { $a->{bin} <=> $b->{bin} } @t_sort;
        printXML($lib, $type, \@range_arr);

        #foreach $bin (@range_arr){
        #    print $bin->{bin}."\t".$bin->{count}."\n";
        #}
    } else {
        print "Nothing to do, empty result set\n";
    }
}

###############################################################################
sub printXML {
    my ($lib, $type, $arr) = @_;

    my $filename = $this->{output_dir} . "/".$type."_HISTOGRAM_".$lib.".xml";

    my $output = new IO::File(">$filename") or die "Could not open file $filename to write\n";
    my $writer = new XML::Writer(OUTPUT=>$output);

    #print $arr[0]."\n";

    $writer->xmlDecl("UTF-8");

    $writer->startTag("root");
    foreach my $bin (@{$arr}){
        $writer->emptyTag("CATEGORY", 'LABEL'=> $bin->{bin}, 'VALUE'=> $bin->{count});
    }
    $writer->endTag("root");
    $writer->end();
    $output->close();
}
