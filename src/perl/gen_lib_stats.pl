#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

  gen_lib_stats.pl

=head1 SYNOPSIS

   USAGE: gen_lib_stats.pl --server server-name --env dbi [--library libraryId]

=head1 OPTIONS

B<--server,-s>
  Server name that need to be processed.

B<--library,-l>
   Specific libraryId whoes taxonomy info to collect

B<--env,-e>
   Specific environment where this script is executed.  Based on these values
   db connection and file locations are set.  Possible values are
   igs, dbi, ageek or test

B<--help,-h>
  This help message


=head1  DESCRIPTION

  This script will process all libraries on a given server.  Get a
  break down of
     ORF types (missing 3', missing 5', incomplete, complete)
     ORF model (bacteria, archaea, phage)
     LIB type (viral only, microbial only, top viral, top microbial)
     FUNCTIONAL and UNClassified

  Counts for each categories are stored in _cnt field, and all sequenceIds
  for each categories are stored in an external file.

=head1  INPUT

  The input is defined with --server which is a domain name only.
     e.g.: calliope (if server name is calliope.dbi.udel.edu)


=head1  OUTPUT

    All counts for each category are stored in "statistics" table on the "server"
    given as input.  All sequenceIds for each category are stored in an
    external file, and its location is stored in db.

=head1  CONTACT

 Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE

 gen_lib_stats.pl --server calliope --env dbi --library 31

=cut

#use lib '/home/zschreib/datastore/docker_testing/opt/ergatis/lib/perl5/';
#use lib '/home/zschreib/perl5/lib/perl5';

use strict;
use warnings;
use DBI;
use Pod::Usage;
use Data::Dumper;
use UTILS_V;
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
if( $options{'help'} ) {
  pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				 'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();
##############################################################################
#### DEFINE GLOBAL VAIRABLES.
##############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

#init db connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{input}", "", "", { RaiseError => 1}) or die $DBI::errstr;

##########################################################################
timer(); #call timer to see when process ended.

print "DEBUG: Setting up sql statements\n";
my $inst_fxn = $dbh->prepare(qq{INSERT INTO statistics (id, libraryId,
							 read_cnt, read_mb,
							 complete_cnt, complete_mb, complete_id,
							 incomplete_cnt, incomplete_mb, incomplete_id,
							 lackstart_cnt, lackstart_mb, lackstart_id,
							 lackstop_cnt, lackstop_mb, lackstop_id,
							 archaea_cnt, archaea_mb, archaea_id,
							 bacteria_cnt, bacteria_mb, bacteria_id,
							 phage_cnt, phage_mb, phage_id,
							 allviral_cnt, allviral_id,
							 topviral_cnt, topviral_id,
							 allmicrobial_cnt, allmicrobial_id,
							 topmicrobial_cnt, topmicrobial_id,
							 fxn_cnt, fxn_id,
							 unassignfxn_cnt, unassignfxn_id,
							 tRNA_cnt, tRNA_id,
							 rRNA_cnt, rRNA_id,
							 orfan_cnt, orfan_id,
							 lineage, deleted, dateCreated)
							 VALUES(?,?,?,?,
							 ?,?,?,
							 ?,?,?,
							 ?,?,?,
							 ?,?,?,
							 ?,?,?,
							 ?,?,?,
							 ?,?,?,
							 ?,?,
							 ?,?,
							 ?,?,
							 ?,?,
							 ?,?,
							 ?,?,
							 ?,?,
							 ?,?,
							 ?,?,
							 ?,?,?)
							});

####get taxonomy
my $tax = $dbh->prepare(qq{SELECT b.domain, count(b.domain)
		       FROM blastp b
		       WHERE b.e_value <= 0.001
				  and b.database_name LIKE 'UNIREF100P'
				  and b.sys_topHit = 1
		       GROUP BY b.domain ORDER BY b.domain desc });

####get all orf sequences
my $all_seq = $dbh->prepare(qq{SELECT distinct s.id,s.header,s.size,s.basepair
			   FROM sequence s
			   WHERE s.typeId=3 });

####get all sig orf hits
my $sig_seq = $dbh->prepare(qq{SELECT distinct b.sequenceId
			  FROM blastp b
			  WHERE b.e_value <= 0.001 });

####get all reads
my $read = $dbh->prepare(qq{SELECT count(s.id), sum(s.size)
			FROM sequence s
			WHERE s.typeId=1});

####get all tRNA
my $tRNA = $dbh->prepare(qq{SELECT t.sequenceId
			FROM tRNA t
			  INNER JOIN sequence s on t.sequenceId = s.id
			WHERE s.typeId=1});

####get all rRNA
my $rRNA = $dbh->prepare(qq{SELECT s.id
			FROM sequence s
			WHERE s.typeId=2});

#### get all top blast hits for a given library
my $top_hits_stmt = qq{SELECT b.sequenceId, MAX(b.db_ranking_code) AS db_ranking_code, b.fxn_topHit
						   FROM blastp b
						   WHERE b.e_value <= 0.001
							and 	(b.database_name = 'UNIREF100P'
								  OR b.database_name = 'METAGENOMES')
							and 	b.sys_topHit=1
                           GROUP BY b.sequenceId
						   ORDER BY b.sequenceId, db_ranking_code desc};

print "DEBUG: Get library in using top_hits_stmt\n";
timer();

my $libId = 1;
my $top_hits_qry = $dbh->prepare($top_hits_stmt);
$top_hits_qry->execute();

########################################
# SPLIT SEQUENCE INTO UNIREF100P ONLY
# AND METAGENOMES ONLY SET AT
# EVALUE CUTOFF AT 0.001
########################################
# get uniref and metagenomes exclusive sequences.
# top_hits_qry returnes only one row per sequence
# db_ranking_code 10=UNIREF100P
#				  10=METAGENOMES
#

#### We can afford to get all data from using using fetchall_hashref('id')
#### because for a given sequence Id we only one record returned
#### if a sequece has a hit to UNIREF record then classify sequence in
#### functional or unclassified functional category, other it will
#### have a metagenome hit.

print STDOUT "DEBUG: Read all top_hits in results hash ref\n";
timer();

my $result = $top_hits_qry->fetchall_hashref('sequenceId');

#### all uniref hits are further divided based on fxnal flag into
#### fxnal_hash or unclassified_hash
my %meta_hash;
my %uniref_hash;
my %fxnal_hash;
my %unclassified_hash;
my %orfan_hash;

print STDOUT "DEBUG: Separeate uniref (fxnal and unclassified) and metagenome hits\n";
timer();

foreach my $k (keys %{$result}) {
	if ($result->{$k}->{db_ranking_code} == 100) {
        $uniref_hash{$k} = $result->{$k};
	} else {
		$meta_hash{$k} = $result->{$k};
	}
}

#### reset result array of hash
$result = undef;

###########################################
# GET FUNCTIONAL AND UNASSIGNED FUNCTIONAL
# CATEGORIES FOR ALL UNIREF100P SEQUENCES
# AT EVALUE CUTOFF OF 0.001
###########################################

# for all uniref only sequences get functional/unassigned protein info.
print STDOUT "check if fxnal hit per uniref record\n";

my @arr = keys %uniref_hash;
my $k = hasFunctionalHit(@arr);

foreach my $sequenceId (keys %uniref_hash) {
    #divide all hits in fxn and unclassified.
    if (exists $k->{$sequenceId}) {
        $fxnal_hash{$sequenceId} = $uniref_hash{$sequenceId}
    }
    else {
        $unclassified_hash{$sequenceId} = $uniref_hash{$sequenceId}
    }
}

print "DEBUG: hasfxnhit " . scalar(keys %{$k})."\n";
print "DEBUG: uniref hash " . scalar(keys %uniref_hash)."\n";
print "DEBUG: fxnal hash " . scalar(keys %fxnal_hash)."\n";
print "DEBUG: unclassified " . scalar(keys %unclassified_hash)."\n";

###############################################
# CALCULATE ENVIRONMENTAL CATEGORIES FOR
# EACH LIBRARY, VIRAL ONLY, TOP VIRAL,
# MICORBIAL ONLY, TOP MICORBIAL
# AT EVALUE CUTOFF AT 0.001
###############################################
my %env=();
$env{'top_viral'}=0;
$env{'top_viral_list'}="";
$env{'top_micro'}=0;
$env{'top_micro_list'}="";
$env{'viral'}=0;
$env{'viral_list'}="";
$env{'micro'}=0;
$env{'micro_list'}="";

print "DEBUG: Gen env info such as top_viral, top_micro etc\n";
timer();

foreach my $seqid (keys %meta_hash) {
	#### get all blast hits for a seq.
	my $sth = $dbh->prepare(qq{SELECT b.id, b.hit_name, b.sys_topHit, b.query_name, b.hit_description
				FROM blastp b
				WHERE b.sequenceId=?
					and b.e_value<=0.001
					and b.database_name='METAGENOMES'
				ORDER BY b.id,b.sys_topHit});
	$sth->execute(($seqid));

	my $top_hit="";
	my $same_hit=1;

	#### loop through all blast results for a sequence
	while (my $row = $sth->fetchrow_hashref) {

		my $lib_type = "micro";
		if ($$row{hit_description} =~ /(viral) metagenome from/i) {
			$lib_type = $1;
		}

		if ($$row{sys_topHit} == 1) {
			if ($lib_type =~ /viral/i) {
				$top_hit = "viral";
			} else { $top_hit = $lib_type; }
		}

		if (($lib_type !~ /$top_hit/i) && ($same_hit)) {
			$same_hit = 0;
		}
	}

	#### viral or microbial only assignment.
	#### if both $top_hit and $other_hits are same or
	#### $other_hit is empyt i.e only one hit
	my $env_type = $top_hit;
	my $env_type_list = $env_type ."_list";

	#### have multiple hits and top hit is different form other hits.
	#### it is possible for $top_hit ne $other_hit if only one hit
	if (!$same_hit) {
		$env_type = "top_" . $top_hit;
		$env_type_list = $env_type ."_list";
	}

	$env{$env_type} += 1;
	if (length($env{$env_type_list})) {
		$env{$env_type_list} .= "," . $seqid;
	} else {
		$env{$env_type_list} = $seqid;
	}
}

 #################################################
 #### CALCULATE TAXONOMY LINEAGE AT DOMAIN LEVEL
 #### AT EVALUE CUTOFF AT 0.001
 #################################################
 #### get domain taxonomy count.
 my ($type,$count,$lineage);
 $tax->execute();
 $tax->bind_col(1,\$type);
 $tax->bind_col(2,\$count);
 $lineage = "";

print "DEBUG: get taxonomy count\n";
timer();

 while($tax->fetch) {
	if (!length($type)) {
		$type = "Unclassified";
	}
	if (length($lineage)) {
		$lineage = $lineage.";".$type.":".$count;
	}
	else { $lineage = $type.":".$count; }
}

 ##################################################
 # CALCULATE ORF CATEGORIES and TYPES
 # FOR EACH LIBRARY AT EVALUE CUTOFF OF 0.001
 ##################################################
 my %orf=();
 $orf{'comp_cnt'}=0;
 $orf{'comp_lst'}="";
 $orf{'comp_mb'}=0;
 $orf{'incomp_cnt'}=0;
 $orf{'incomp_lst'}="";
 $orf{'incomp_mb'}=0;
 $orf{'start_cnt'}=0;
 $orf{'start_lst'}="";
 $orf{'start_mb'}=0;
 $orf{'stop_cnt'}=0;
 $orf{'stop_lst'}="";
 $orf{'stop_mb'}=0;
 $orf{'bacteria_cnt'}=0;
 $orf{'bacteria_lst'}="";
 $orf{'bacteria_mb'}=0;
 $orf{'archaea_cnt'}=0;
 $orf{'archaea_lst'}="";
 $orf{'archaea_mb'}=0;
 $orf{'phage_cnt'}=0;
 $orf{'phage_lst'}="";
 $orf{'phage_mb'}=0;

print "DEBUG: get complete, incomplete and other orf types including orfans\n";
timer();

$all_seq->execute();
my $all_seq_hash = $all_seq->fetchall_hashref('id');

foreach my $k (keys %{$all_seq_hash}) {
    my %opts;

	map { $opts{$1} = $2 if( /([^=]+)\s*=\s*([^=]+)/ ) } split(/\s+/, $all_seq_hash->{$k}->{header});

	if ($opts{type} =~ /lack[_|\s]stop/i){
		$opts{type} = "stop";
	} elsif ($opts{type} =~ /lack[_|\s]start/i){
		$opts{type} = "start";
	} elsif ($opts{type} =~ /incomplete/i){
		$opts{type} = "incomp";
	} elsif ($opts{type} =~ /complete/i){
		$opts{type} = "comp";
	}

	#### set model stats.
	$orf{$opts{model}.'_cnt'}++;

	#if * at the end of bases, don't count it
	if ($all_seq_hash->{$k}->{basepair} =~ /\*$/) {
		$orf{$opts{model}.'_mb'} += ($all_seq_hash->{$k}->{size}-1);
	} else {
		$orf{$opts{model}.'_mb'} += $all_seq_hash->{$k}->{size};
	}

	if (length($orf{$opts{model}.'_lst'})) {
		$orf{$opts{model}.'_lst'} = $orf{$opts{model}.'_lst'} . "," . $all_seq_hash->{$k}->{id};
	} else { $orf{$opts{model}.'_lst'} = $all_seq_hash->{$k}->{id}; }

	#### set type stats.
	$orf{$opts{type}.'_cnt'}++;

	####do not count * in bases
	if ($all_seq_hash->{$k}->{basepair} =~ /\*$/) {
		$orf{$opts{type}.'_mb'} += ($all_seq_hash->{$k}->{size} - 1);
	} else { $orf{$opts{type}.'_mb'} += $all_seq_hash->{$k}->{size}; }

	if (length($orf{$opts{type}.'_lst'})) {
		$orf{$opts{type}.'_lst'} = $orf{$opts{type}.'_lst'} . "," . $all_seq_hash->{$k}->{id};
	} else { $orf{$opts{type}.'_lst'} = $all_seq_hash->{$k}->{id}; }

    #### generate orfan list by comparing all uniref/metagenome sequences that have
    #### a blast hit that is significant that e_value < 0.001 if no such seq exist in
    #### blastp table then classify the sequence as orfan.

    #### if sequence id does not exists in meta_hash, fxnal_hash or unclassified_hash
    #### then seq is an orfan
    unless ((exists $meta_hash{$k}) or (exists $fxnal_hash{$k}) or (exists $unclassified_hash{$k})) {
        $orfan_hash{$k} = $all_seq_hash->{$k};
    }
}

##################################################
#### GET READ COUNT AND MEGABASES
##################################################
print "DEBUG: get read counts\n";
timer();

$read->execute();
my $read_row = $read->fetchall_arrayref([],1);
my %read_s = ();
$read_s{'count'} = $read_row->[0]->[0];
$read_s{'mb'} = $read_row->[0]->[1];

##################################################
#### GET TRNA COUNT AND LIST
##################################################
print "DEBUG: get trna\n";
timer();

$tRNA->execute();
my %tRNA_s = ();
$tRNA_s{'count'} = 0;
$tRNA_s{'lst'} = "";

while (my @rslt = $tRNA->fetchrow_array()) {
	$tRNA_s{'count'}++;
	if (length($tRNA_s{'lst'})) {
		$tRNA_s{'lst'} = $tRNA_s{'lst'} . "," . $rslt[0];
	} else { $tRNA_s{'lst'} = $rslt[0]; }
}

###################################################
#### GET RRNA COUNT AND LIST
###################################################
print "DEBUG: get rrna\n";
timer();

$rRNA->execute();
my %rRNA_s = ();
$rRNA_s{'count'} = 0;
$rRNA_s{'lst'} = "";

while (my @rslt = $rRNA->fetchrow_array()) {
	$rRNA_s{'count'}++;
	if (length($rRNA_s{'lst'})) {
		$rRNA_s{'lst'} = $rRNA_s{'lst'} . "," . $rslt[0];
	} else { $rRNA_s{'lst'} = $rslt[0]; }
}


#################################################
#### INSERT STATISTICS INTO TABLE
#################################################

####create output file for functional_id list and unclass_id list
my $output_dir = "$options{outdir}/idFiles";

my %file_output_list= (
    "functional_list" => $output_dir."/fxnIdList_".$libId.".txt",
    "unclassified_list" => $output_dir."/unClassIdList_".$libId.".txt",
    "viral_list" => $output_dir."/viralList_".$libId.".txt",
    "top_viral_list" => $output_dir."/topViralList_".$libId.".txt",
    "micro_list" => $output_dir."/microList_".$libId.".txt",
    "top_micro_list" => $output_dir."/topMicroList_".$libId.".txt",
    "comp_lst" => $output_dir."/compORFList_".$libId.".txt",
    "incomp_lst" => $output_dir."/incompORFList_".$libId.".txt",
    "start_lst" => $output_dir."/startORFList_".$libId.".txt",
    "stop_lst" => $output_dir."/stopORFList_".$libId.".txt",
    "archaea_lst" => $output_dir."/arcORFList_".$libId.".txt",
    "bacteria_lst" => $output_dir."/bacORFList_".$libId.".txt",
    "phage_lst" => $output_dir."/phgORFList_".$libId.".txt",
    "orfan" => $output_dir."/orfanList_".$libId.".txt",
    "tRNA" => $output_dir."/tRNAList_".$libId.".txt",
    "rRNA" => $output_dir."/rRNAList_".$libId.".txt"
);

print "DEBUG: Gather data to insert into table\n";
timer();

foreach my $key (keys %file_output_list) {

   open(OUT, ">", $file_output_list{$key} ) or
	  $logger->logdie("Could not open file $file_output_list{$key} to write");

   if ($key =~ /functional_list/i) {
	  print OUT join(",", keys %fxnal_hash);
   } elsif ($key =~ /unclassified_list/i) {
	  print OUT join(",", keys %unclassified_hash);
   } elsif ($key =~ /viral_list|top_viral_list|micro_list|top_micro_list/i) {
	  print OUT $env{$key};
   } elsif ($key =~ /comp_lst|incomp_lst|start_lst|stop_lst|archaea_lst|bacteria_lst|phage_lst/i) {
	  print OUT $orf{$key};
   } elsif ($key =~ /orfan/i){
	  print OUT join(",", keys %orfan_hash);
   } elsif ($key =~ /tRNA/i){
	  print OUT $tRNA_s{'lst'};
   } elsif ($key =~ /rRNA/i){
	  print OUT $rRNA_s{'lst'};
   }

   close OUT;
}

print "DEBUG: get max row id\n";
timer();

my $max_id = 0;
my $sth = $dbh->prepare('SELECT max(id) as max_id FROM statistics');
$sth->execute;

while ( my $result = $sth->fetchrow_hashref() ) {
	if (defined $result->{max_id}) {
		$max_id = $result->{max_id} + 1;
	} else {
		$max_id = 1;
	}
}

print "DEBUG: run sql insert\n";
timer();

$inst_fxn->execute($max_id, $libId, $read_s{count}, $read_s{mb},
				   $orf{comp_cnt}, $orf{comp_mb}, "compORFList_${libId}.txt",
				   $orf{incomp_cnt}, $orf{incomp_mb}, "incompORFList_${libId}.txt",
				   $orf{start_cnt}, $orf{start_mb}, "startORFList_${libId}.txt",
				   $orf{stop_cnt}, $orf{stop_mb}, "stopORFList_${libId}.txt",
				   $orf{archaea_cnt}, $orf{archaea_mb}, "arcORFList_${libId}.txt",
				   $orf{bacteria_cnt}, $orf{bacteria_mb}, "bacORFList_${libId}.txt",
				   $orf{phage_cnt}, $orf{phage_mb}, "phgORFList_${libId}.txt",
				   $env{viral}, "viralList_${libId}.txt",
				   $env{top_viral}, "topViralList_${libId}.txt",
				   $env{micro}, "microList_${libId}.txt",
				   $env{top_micro}, "topMicroList_${libId}.txt",
				   scalar(keys %fxnal_hash), "fxnIdList_${libId}.txt",
				   scalar(keys %unclassified_hash), "unClassIdList_${libId}.txt",
				   $tRNA_s{count}, "tRNAList_${libId}.txt",
				   $rRNA_s{count}, "rRNAList_${libId}.txt",
				   scalar(keys %orfan_hash), "orfanList_${libId}.txt",
				   $lineage, 0, getTimeStamp()
				   );
$max_id++;

print "DEBUG: Complete\n";
timer();

$logger->info("General library stats completed");
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
sub hasFunctionalHit {
   my $seqId = join(",", @_);

   my $qry_stmt = qq{SELECT distinct b.sequenceId, b.fxn_topHit
								  FROM blastp b
								  WHERE b.fxn_topHit=1
									and b.e_value<=0.001
									and b.sequenceId in ($seqId)};

   my $fxn_hit = $dbh->prepare($qry_stmt);
   $fxn_hit->execute();

   my $result = $fxn_hit->fetchall_hashref('sequenceId');

   return $result;
}

###############################################################################
sub getTimeStamp {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d-%02d-%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
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
