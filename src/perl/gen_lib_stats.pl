#! /usr/bin/perl

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
my $results =
  GetOptions( \%options, 'input|i=s', 'outdir|o=s', 'log|l=s', 'debug|d=s',
    'help|h' )
  || pod2usage();

## display documentation
if ( $options{'help'} ) {
    pod2usage( { -exitval => 0, -verbose => 2, -output => \*STDERR } );
}

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger(
    'LOG_FILE'  => $logfile,
    'LOG_LEVEL' => $options{'debug'}
);
$logger = $logger->get_logger();
##############################################################################
#### DEFINE GLOBAL VAIRABLES.
##############################################################################
## make sure everything passed was peachy
&check_parameters( \%options );

#init db connection
my $dbh = DBI->connect( "dbi:SQLite:dbname=$options{input}",
    "", "", { RaiseError => 1 } )
  or die $DBI::errstr;

##########################################################################
timer();    #call timer to see when process ended.

print STDOUT "Setup all SQL queries\n";

my $inst_fxn = $dbh->prepare(
    qq{INSERT INTO statistics (id, libraryId,
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
							}
);

####get taxonomy
my $tax = $dbh->prepare(
    qq{SELECT b.domain, count(b.domain)
		       FROM blastp b
				INNER JOIN sequence s on s.id = b.sequenceId
		       WHERE s.libraryId = ?
				  and b.e_value <= 0.001
				  and b.database_name LIKE 'UNIREF100P'
				  and b.sys_topHit = 1
		       GROUP BY b.domain ORDER BY b.domain desc }
);

####get all orf sequences
my $all_seq = $dbh->prepare(
    qq{SELECT distinct s.id,s.header,s.size,s.basepair
			   FROM sequence s
				INNER JOIN sequence_relationship sr on s.id = sr.objectId
			   WHERE s.libraryId=?
			     and sr.typeId=3 }
);

####get all sig orf hits
my $sig_seq = $dbh->prepare(
    qq{SELECT distinct s.id, s.header, s.size
			  FROM blastp b
				INNER JOIN sequence s on b.sequenceId = s.id
				INNER JOIN sequence_relationship sr on s.id = sr.objectId
			  WHERE s.libraryId=?
			    and b.e_value <= 0.001
			    and sr.typeId = 3
			    and (b.database_name like 'UNIREF100P' OR
					 b.database_name like 'METAGENOMES') }
);

####get all reads
my $read = $dbh->prepare(
    qq{SELECT count(s.id), sum(s.size)
			FROM sequence s
			  INNER JOIN sequence_relationship sr on s.id = sr.objectId
			WHERE s.libraryId=?
			  and sr.typeId=1}
);

####get all tRNA
my $tRNA = $dbh->prepare(
    qq{SELECT t.sequenceId
			FROM tRNA t
			  INNER JOIN sequence s on t.sequenceId = s.id
			  INNER JOIN sequence_relationship sr on s.id = sr.objectId
			WHERE s.libraryId=?
			  and sr.typeId=1}
);

####get all rRNA
my $rRNA = $dbh->prepare(
    qq{SELECT s.id
			FROM sequence s
			  INNER JOIN sequence_relationship sr on s.id = sr.objectId
			WHERE s.libraryId=?
			  and sr.typeId=2}
);

#### get all top blast hits for a given library
my $top_hits_stmt =
  qq{SELECT b.sequenceId, MAX(b.db_ranking_code) AS db_ranking_code
                        FROM blastp b
                        WHERE b.e_value <= 0.001
                         and 	(b.database_name = 'UNIREF100P'
                               OR b.database_name = 'METAGENOMES')
                         and 	b.sys_topHit=1
                        GROUP BY b.sequenceId
                        ORDER BY b.sequenceId, db_ranking_code desc};

                        # qq{SELECT b.sequenceId, MAX(b.db_ranking_code) AS db_ranking_code
                      	# 					   FROM blastp b INNER JOIN sequence s ON s.id=b.sequenceId
                      	# 					   WHERE s.libraryId = ?
                      	# 						and 	b.e_value <= 0.001
                      	# 						and 	(b.database_name = 'UNIREF100P'
                      	# 							  OR b.database_name = 'METAGENOMES')
                      	# 						and 	b.sys_topHit=1
                      	# 					   GROUP BY b.sequenceId
                      	# 					   ORDER BY b.sequenceId, db_ranking_code desc};

print STDOUT "Start GenLibStats\n";

print STDOUT "Read sequence table\n";
my $libId        = 1;
my $top_hits_qry = $dbh->prepare($top_hits_stmt);
$top_hits_qry->execute($libId);

print STDOUT "Read all blast tophit results for UNIREF and MGOL\n";

#### full all records for top_hits_qry instead of prcessing one row at a time
#### could help with performance as all data will be read in memory,
#### and processed in memory.
my $result = $top_hits_qry->fetchall_arrayref( {} );

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
print STDOUT "Separeate uniref and metagenome hits\n";

my ( @uniref_arr, @meta_arr );
my $prev = -1;

#### commen of per row process to test if all rows in memory
#### work faster than this
# while (my $result = $top_hits_qry->fetchrow_hashref()) {
# 	if ($$result{db_ranking_code} == 100) {
# 		push @uniref_arr, $$result{sequenceId};
# 	} else {
# 		push @meta_arr , $$result{sequenceId};
# 	}
# }

foreach my $row ( @{$result} ) {
    if ( $row->{db_ranking_code} == 100 ) {
        push @uniref_arr, $row->{sequenceId};
    }
    else {
        push @meta_arr, $row->{sequenceId};
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
my $functional_count   = 0;
my $functional_list    = "";
my $unclassified_count = 0;
my $unclassified_list  = "";
my $null_str           = "NULL";

print STDOUT "check if fxnal hit per uniref record\n";

foreach my $sequenceId (@uniref_arr) {

    #divide all hits in fxn and unclassified.
    if ( hasFunctionalHit($sequenceId) ) {
        $functional_count++;
        $functional_list .= $sequenceId . ",";
    }
    else {
        $unclassified_count++;
        $unclassified_list .= $sequenceId . ",";
    }
}

#remove last comma
$functional_list =~ s/,$//;
$unclassified_list =~ s/,$//;

###############################################
# CALCULATE ENVIRONMENTAL CATEGORIES FOR
# EACH LIBRARY, VIRAL ONLY, TOP VIRAL,
# MICORBIAL ONLY, TOP MICORBIAL
# AT EVALUE CUTOFF AT 0.001
###############################################
my %env = ();
$env{'top_viral'}      = 0;
$env{'top_viral_list'} = "";
$env{'top_micro'}      = 0;
$env{'top_micro_list'} = "";
$env{'viral'}          = 0;
$env{'viral_list'}     = "";
$env{'micro'}          = 0;
$env{'micro_list'}     = "";

print STDOUT "Get all viral, microbial list\n";

foreach my $seqid (@meta_arr) {
    #### get all blast hits for a seq.
    my $sth = $dbh->prepare(
            qq{SELECT b.id, b.hit_name, b.sys_topHit, b.query_name, b.hit_description
				FROM blastp b
				WHERE b.sequenceId=?
					and b.e_value<=0.001
					and b.database_name='METAGENOMES'
				ORDER BY b.id,b.sys_topHit}
    );
    $sth->execute( ($seqid) );
    $result = $sth->fetchall_arrayref( {} );

    my $top_hit  = "";
    my $same_hit = 1;

    #### loop through all blast results for a sequence
    # while (my $row = $sth->fetchrow_hashref) {

    #### try fetchall_arrayref instead of fetchrow and see if all
    #### in memory process will increase the speed performance
    foreach my $row ( @{$result} ) {
        my $lib_type = "micro";
        if ( $row->{hit_description} =~ /(viral) metagenome from/i ) {
            $lib_type = $1;
        }

        if ( $row->{sys_topHit} == 1 ) {
            if ( $lib_type =~ /viral/i ) {
                $top_hit = "viral";
            }
            else {
                $top_hit = $lib_type;
            }
        }

        if ( ( $lib_type !~ /$top_hit/i ) && ($same_hit) ) {
            $same_hit = 0;
        }
    }

    #### reset result to release memory requirement.
    $result = undef;

    #### viral or microbial only assignment.
    #### if both $top_hit and $other_hits are same or
    #### $other_hit is empyt i.e only one hit
    my $env_type      = $top_hit;
    my $env_type_list = $env_type . "_list";

    #### have multiple hits and top hit is different form other hits.
    #### it is possible for $top_hit ne $other_hit if only one hit
    if ( !$same_hit ) {
        $env_type      = "top_" . $top_hit;
        $env_type_list = $env_type . "_list";
    }

    $env{$env_type} += 1;
    if ( length( $env{$env_type_list} ) ) {
        $env{$env_type_list} .= "," . $seqid;
    }
    else {
        $env{$env_type_list} = $seqid;
    }
}

#################################################
#### CALCULATE TAXONOMY LINEAGE AT DOMAIN LEVEL
#### AT EVALUE CUTOFF AT 0.001
#################################################
#### get domain taxonomy count.
my ( $type, $count, $lineage );
$tax->execute( ($libId) );
$tax->bind_col( 1, \$type );
$tax->bind_col( 2, \$count );
$lineage = "";

print STDOUT "Get Taxonomy lineage\n";

while ( $tax->fetch ) {
    if ( !length($type) ) {
        $type = "Unclassified";
    }
    if ( length($lineage) ) {
        $lineage = $lineage . ";" . $type . ":" . $count;
    }
    else { $lineage = $type . ":" . $count; }
}

##################################################
# CALCULATE ORF CATEGORIES and TYPES
# FOR EACH LIBRARY AT EVALUE CUTOFF OF 0.001
##################################################
my %orf = ();
$orf{'comp_cnt'}     = 0;
$orf{'comp_lst'}     = "";
$orf{'comp_mb'}      = 0;
$orf{'incomp_cnt'}   = 0;
$orf{'incomp_lst'}   = "";
$orf{'incomp_mb'}    = 0;
$orf{'start_cnt'}    = 0;
$orf{'start_lst'}    = "";
$orf{'start_mb'}     = 0;
$orf{'stop_cnt'}     = 0;
$orf{'stop_lst'}     = "";
$orf{'stop_mb'}      = 0;
$orf{'bacteria_cnt'} = 0;
$orf{'bacteria_lst'} = "";
$orf{'bacteria_mb'}  = 0;
$orf{'archaea_cnt'}  = 0;
$orf{'archaea_lst'}  = "";
$orf{'archaea_mb'}   = 0;
$orf{'phage_cnt'}    = 0;
$orf{'phage_lst'}    = "";
$orf{'phage_mb'}     = 0;

print STDOUT "Get ORF count and ORF type list\n";

$all_seq->execute( ($libId) );
$result = $all_seq->fetchall_arrayref( {} );

#  while (my $row = $all_seq->fetchrow_hashref) {
foreach my $row ( @{$result} ) {
    my %opts;

    map { $opts{$1} = $2 if (/([^=]+)\s*=\s*([^=]+)/) }
      split( /\s+/, $row->{header} );

    if ( $opts{type} =~ /lack[_|\s]stop/i ) {
        $opts{type} = "stop";
    }
    elsif ( $opts{type} =~ /lack[_|\s]start/i ) {
        $opts{type} = "start";
    }
    elsif ( $opts{type} =~ /incomplete/i ) {
        $opts{type} = "incomp";
    }
    elsif ( $opts{type} =~ /complete/i ) {
        $opts{type} = "comp";
    }

    #### set model stats.
    $orf{ $opts{model} . '_cnt' }++;

    #if * at the end of bases, don't count it
    if ( $row->{basepair} =~ /\*$/ ) {
        $orf{ $opts{model} . '_mb' } += ( $row->{size} - 1 );
    }
    else {
        $orf{ $opts{model} . '_mb' } += $row->{size};
    }

    if ( length( $orf{ $opts{model} . '_lst' } ) ) {
        $orf{ $opts{model} . '_lst' } =
          $orf{ $opts{model} . '_lst' } . "," . $row->{id};
    }
    else {
        $orf{ $opts{model} . '_lst' } = $row->{id};
    }

    #### set type stats.
    $orf{ $opts{type} . '_cnt' }++;

    ####do not count * in bases
    if ( $row->{basepair} =~ /\*$/ ) {
        $orf{ $opts{type} . '_mb' } += ( $row->{size} - 1 );
    }
    else {
        $orf{ $opts{type} . '_mb' } += $row->{size};
    }

    if ( length( $orf{ $opts{type} . '_lst' } ) ) {
        $orf{ $opts{type} . '_lst' } =
          $orf{ $opts{type} . '_lst' } . "," . $row->{id};
    }
    else {
        $orf{ $opts{type} . '_lst' } = $row->{id};
    }
}

#### reset result var to release memory requirement
$result = undef;

##################################################
#### GET ORFAN COUNT AT EVALUE CUTOFF OF 0.001
##################################################
my $sigcnt = 0;
my $siglst = "";

print STDOUT "Get ORFan count\n";

$sig_seq->execute( ($libId) );
$result = $sig_seq->fetchall_arrayref( {} );

foreach my $val ( @{$result} ) {
    $sigcnt++;
    #### get all significant hit ids.
    if ( length($siglst) ) {
        $siglst .= "," . $val->{id};
    }
    else { $siglst = $val->{id}; }
}

#### reset result variable
$result = undef;

$all_seq->execute( ($libId) );
$result = $all_seq->fetchall_arrayref( {} );

my $allcount = 0;
my %orfan    = ();

$orfan{'count'} = 0;
$orfan{'lst'}   = "";

foreach my $val ( @{$result} ) {
    $allcount++;
    if ( index( $siglst, $val->{id} ) < 0 ) {
        $orfan{'count'}++;

        if ( length( $orfan{'lst'} ) ) {
            $orfan{'lst'} .= "," . $val->{id};
        }
        else { $orfan{'lst'} = $val->{id}; }
    }
}

##################################################
#### GET READ COUNT AND MEGABASES
##################################################
$read->execute( ($libId) );

print STDOUT "Get read count\n";

my $read_row = $read->fetchall_arrayref( [], 1 );
my %read_s = ();

$read_s{'count'} = $read_row->[0]->[0];
$read_s{'mb'}    = $read_row->[0]->[1];

##################################################
#### GET TRNA COUNT AND LIST
##################################################
$tRNA->execute( ($libId) );
my %tRNA_s = ();
$tRNA_s{'count'} = 0;
$tRNA_s{'lst'}   = "";

print STDOUT "Get tRNA count\n";

$result = $tRNA->fetchall_arrayref({});

foreach my $rslt ( @{$result} ) {
    if ( length( $tRNA_s{'lst'} ) ) {
        $tRNA_s{'lst'} .= "," . $rslt->{sequenceId};
    }
    else {
        $tRNA_s{'lst'} = $rslt->{sequenceId};
    }
    $tRNA_s{'count'}++;
}

#### reset results variable
$result = undef;

###################################################
#### GET RRNA COUNT AND LIST
###################################################
$rRNA->execute( ($libId) );
my %rRNA_s = ();
$rRNA_s{'count'} = 0;
$rRNA_s{'lst'}   = "";

print STDOUT "Get rRNA count\n";

while ( my @rslt = $rRNA->fetchrow_array() ) {
    $rRNA_s{'count'}++;
    if ( length( $rRNA_s{'lst'} ) ) {
        $rRNA_s{'lst'} = $rRNA_s{'lst'} . "," . $rslt[0];
    }
    else { $rRNA_s{'lst'} = $rslt[0]; }
}

#################################################
#### INSERT STATISTICS INTO TABLE
#################################################

####create output file for functional_id list and unclass_id list
my $output_dir = "$options{outdir}/idFiles";

my %file_output_list = (
    "functional_list"   => $output_dir . "/fxnIdList_" . $libId . ".txt",
    "unclassified_list" => $output_dir . "/unClassIdList_" . $libId . ".txt",
    "viral_list"        => $output_dir . "/viralList_" . $libId . ".txt",
    "top_viral_list"    => $output_dir . "/topViralList_" . $libId . ".txt",
    "micro_list"        => $output_dir . "/microList_" . $libId . ".txt",
    "top_micro_list"    => $output_dir . "/topMicroList_" . $libId . ".txt",
    "comp_lst"          => $output_dir . "/compORFList_" . $libId . ".txt",
    "incomp_lst"        => $output_dir . "/incompORFList_" . $libId . ".txt",
    "start_lst"         => $output_dir . "/startORFList_" . $libId . ".txt",
    "stop_lst"          => $output_dir . "/stopORFList_" . $libId . ".txt",
    "archaea_lst"       => $output_dir . "/arcORFList_" . $libId . ".txt",
    "bacteria_lst"      => $output_dir . "/bacORFList_" . $libId . ".txt",
    "phage_lst"         => $output_dir . "/phgORFList_" . $libId . ".txt",
    "orfan"             => $output_dir . "/orfanList_" . $libId . ".txt",
    "tRNA"              => $output_dir . "/tRNAList_" . $libId . ".txt",
    "rRNA"              => $output_dir . "/rRNAList_" . $libId . ".txt"
);

print STDOUT "Insert data into stats table\n";

foreach my $key ( keys %file_output_list ) {

    open( OUT, ">", $file_output_list{$key} )
      or
      $logger->logdie("Could not open file $file_output_list{$key} to write");

    if ( $key =~ /functional_list/i ) {
        print OUT $functional_list;
    }
    elsif ( $key =~ /unclassified_list/i ) {
        print OUT $unclassified_list;
    }
    elsif ( $key =~ /viral_list|top_viral_list|micro_list|top_micro_list/i ) {
        print OUT $env{$key};
    }
    elsif ( $key =~/comp_lst|incomp_lst|start_lst|stop_lst|archaea_lst|bacteria_lst|phage_lst/i )
    {
        print OUT $orf{$key};
    }
    elsif ( $key =~ /orfan/i ) {
        print OUT $orfan{'lst'};
    }
    elsif ( $key =~ /tRNA/i ) {
        print OUT $tRNA_s{'lst'};
    }
    elsif ( $key =~ /rRNA/i ) {
        print OUT $rRNA_s{'lst'};
    }

    close OUT;
}

my $max_id = 0;
my $sth    = $dbh->prepare('SELECT max(id) as max_id FROM statistics');
$sth->execute;

while ( my $result = $sth->fetchrow_hashref() ) {
    if ( defined $result->{max_id} ) {
        $max_id = $result->{max_id} + 1;
    }
    else {
        $max_id = 1;
    }
}

$inst_fxn->execute(
    $max_id,                      $libId,
    $read_s{count},               $read_s{mb},
    $orf{comp_cnt},               $orf{comp_mb},
    "compORFList_${libId}.txt",   $orf{incomp_cnt},
    $orf{incomp_mb},              "incompORFList_${libId}.txt",
    $orf{start_cnt},              $orf{start_mb},
    "startORFList_${libId}.txt",  $orf{stop_cnt},
    $orf{stop_mb},                "stopORFList_${libId}.txt",
    $orf{archaea_cnt},            $orf{archaea_mb},
    "arcORFList_${libId}.txt",    $orf{bacteria_cnt},
    $orf{bacteria_mb},            "bacORFList_${libId}.txt",
    $orf{phage_cnt},              $orf{phage_mb},
    "phgORFList_${libId}.txt",    $env{viral},
    "viralList_${libId}.txt",     $env{top_viral},
    "topViralList_${libId}.txt",  $env{micro},
    "microList_${libId}.txt",     $env{top_micro},
    "topMicroList_${libId}.txt",  $functional_count,
    "fxnIdList_${libId}.txt",     $unclassified_count,
    "unClassIdList_${libId}.txt", $tRNA_s{count},
    "tRNAList_${libId}.txt",      $rRNA_s{count},
    "rRNAList_${libId}.txt",      $orfan{count},
    "orfanList_${libId}.txt",     $lineage,
    0,                            getTimeStamp()
);
$max_id++;

$logger->info("General library stats for $options{sample} completed");
exit(0);

###############################################################################
####  SUBS
###############################################################################

sub check_parameters {
    my $options = shift;

    my @required = qw(input outdir);

    foreach my $key (@required) {
        unless ( $options{$key} ) {
            pod2usage(
                {
                    -exitval => 2,
                    -message => "ERROR: Input $key not defined",
                    -verbose => 1,
                    -output  => \*STDERR
                }
            );
            $logger->logdie("No input defined, plesae read perldoc $0\n\n");
        }
    }

    make_path( $options{outdir} . "/idFiles" );
    make_path( $options{outdir} . "/xDocs" );
}

###############################################################################
sub hasFunctionalHit {
    my $seqId = shift;

    my $fxn_hit = $dbh->prepare(
        qq{SELECT	b.id
								  FROM		blastp b
								  WHERE 	b.fxn_topHit=1
									and 	b.e_value<=0.001
									and		b.sequenceId=?}
    );
    $fxn_hit->execute($seqId);

    while ( my $hits = $fxn_hit->fetchrow_hashref() ) {
        return 1;
    }

    return 0;
}

###############################################################################
sub getTimeStamp {

    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
      localtime(time);
    my $nice_timestamp = sprintf(
        "%04d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $mday, $hour, $min, $sec
    );
    return $nice_timestamp;
}

###############################################################################
sub timer {
    my @months   = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my (
        $second,     $minute,    $hour,
        $dayOfMonth, $month,     $yearOffset,
        $dayOfWeek,  $dayOfYear, $daylightSavings
    ) = localtime();
    my $year = 1900 + $yearOffset;
    my $theTime =
"$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
    print "Time now: " . $theTime . "\n";
}
