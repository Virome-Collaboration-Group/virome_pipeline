#!/usr/bin/perl

=head1 NAME
   init_db.pl

=head1 SYNOPSIS

    USAGE: init-db.pl -o=/path/to/output/dir

=head1 OPTIONS

B<--outdir,-o>
   Absolute path to the output repository. (REQUIRED)

B<--log,-l>
   Location of the log file. (Optional)

B<--debug,-d>
   Debugging level of the log file:
     0 = OFF,
     1 = FATAL,
     2 = ERROR,
     3 = WARN,
     4 = INFO,
     5 = DEBUG,
     6 = ALL

B<--help,-h>
   This help message (Optional)

=head1  DESCRIPTION
    Initializes the SQLite database and tables.

=head1  INPUT


=head1  OUTPUT
    SQLite DB file in the output_repository

=head1  CONTACT


==head1 EXAMPLE
   init_db.pl -o=/path/to/output/dir

=cut

use strict;
use warnings;
use DBI;
use DBD::SQLite;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
use UTILS_V;

BEGIN {
    use Ergatis::Logger;
}

##############################################################################
my %options = ();
my $results = GetOptions (\%options,
						  'outdir|o=s',
                          'pstore|p=s',
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

###############################################################################

my $storage = `df -h`;
print STDERR $storage;
print $storage;

## Init the SQLite database file
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{pstore}/processing.sqlite3", "", "", { RaiseError => 1}) or die $DBI::errstr;

$dbh->do( "CREATE TABLE `blastp` (
  `id` INTEGER PRIMARY KEY AUTOINCREMENT,
  `sequenceId` bigint(19) NOT NULL,
  `query_name` varchar(255) DEFAULT NULL,
  `query_length` mediumint(8) NOT NULL,
  `algorithm` varchar(12) NOT NULL,
  `database_name` varchar(15) NOT NULL,
  `db_ranking_code` int(2) NOT NULL DEFAULT '0',
  `hit_name` text,
  `hit_description` text,
  `qry_start` mediumint(8) NOT NULL DEFAULT '0',
  `qry_end` mediumint(8) NOT NULL DEFAULT '0',
  `hit_start` mediumint(8) NOT NULL DEFAULT '0',
  `hit_end` mediumint(8) NOT NULL DEFAULT '0',
  `percent_identity` double NOT NULL DEFAULT '0',
  `percent_similarity` double NOT NULL DEFAULT '0',
  `raw_score` smallint(5) NOT NULL DEFAULT '0',
  `bit_score` double NOT NULL DEFAULT '0',
  `blast_frame` tinyint(1) NOT NULL DEFAULT '0',
  `qry_strand` varchar(10) DEFAULT NULL,
  `subject_length` mediumint(8) NOT NULL DEFAULT '0',
  `e_value` double NOT NULL,
  `domain` varchar(15) DEFAULT NULL,
  `kingdom` varchar(75) DEFAULT NULL,
  `phylum` varchar(75) DEFAULT NULL,
  `class` varchar(75) DEFAULT NULL,
  `order` varchar(75) DEFAULT NULL,
  `family` varchar(75) DEFAULT NULL,
  `genus` varchar(75) DEFAULT NULL,
  `species` varchar(75) DEFAULT NULL,
  `organism` varchar(125) DEFAULT NULL,
  `sys_topHit` tinyint(1) NOT NULL DEFAULT '0',
  `fxn_topHit` tinyint(1) NOT NULL DEFAULT '0',
  `dateCreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `dateModified` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `deleted` tinyint(1) NOT NULL DEFAULT '0'
  );");

$dbh->do( "CREATE INDEX blastp_query_name ON blastp(`query_name`);" );
$dbh->do( "CREATE INDEX blastp_e_value ON blastp(`e_value`);" );
$dbh->do( "CREATE INDEX blastp_sequenceId ON blastp(`sequenceId`);" );
$dbh->do( "CREATE INDEX blastp_databasename ON blastp(`database_name`);" );
$dbh->do( "CREATE INDEX blastp_sysTopHit ON blastp(`sys_topHit`);" );
$dbh->do( "CREATE INDEX blastp_fxnTopHit ON blastp(`fxn_topHit`);" );
$dbh->do( "CREATE INDEX blastp_hitName ON blastp(`hit_name`);" );
#$dbh->do( "CREATE INDEX blastp_query_name ON blastp(query_name,e_value);" );
#$dbh->do( "CREATE INDEX blastp_domain ON blastp(`domain`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX blastp_kingdom ON blastp(`kingdom`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX blastp_phylum ON blastp(`phylum`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX blastp_class ON blastp(`class`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX `blastp_order` ON blastp(`order`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX blastp_family ON blastp(`family`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX blastp_genus ON blastp(`genus`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX blastp_species ON blastp(`species`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX blastp_organism ON blastp(`organism`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX blastp_hit_description ON blastp(`hit_description`,`deleted`,`e_value`);" );
#$dbh->do( "CREATE INDEX blastp_qn_hn ON blastp(`query_name`,`hit_name`,`qry_start`,`hit_start`);" );

## create table blastn
$dbh->do( "CREATE TABLE `blastn` (
  `id` INTEGER PRIMARY KEY AUTOINCREMENT,
  `sequenceId` bigint(19) NOT NULL,
  `query_name` varchar(255) DEFAULT NULL,
  `query_length` mediumint(8) NOT NULL,
  `algorithm` varchar(12) NOT NULL,
  `database_name` varchar(15) NOT NULL,
  `db_ranking_code` int(2) NOT NULL DEFAULT '0',
  `hit_name` text,
  `hit_description` text,
  `qry_start` mediumint(8) NOT NULL DEFAULT '0',
  `qry_end` mediumint(8) NOT NULL DEFAULT '0',
  `hit_start` mediumint(8) NOT NULL DEFAULT '0',
  `hit_end` mediumint(8) NOT NULL DEFAULT '0',
  `percent_identity` double NOT NULL DEFAULT '0',
  `percent_similarity` double NOT NULL DEFAULT '0',
  `raw_score` smallint(5) NOT NULL DEFAULT '0',
  `bit_score` double NOT NULL DEFAULT '0',
  `blast_frame` tinyint(1) NOT NULL DEFAULT '0',
  `qry_strand` varchar(10) DEFAULT NULL,
  `subject_length` mediumint(8) NOT NULL DEFAULT '0',
  `e_value` double NOT NULL,
  `domain` varchar(15) DEFAULT NULL,
  `kingdom` varchar(75) DEFAULT NULL,
  `phylum` varchar(75) DEFAULT NULL,
  `class` varchar(75) DEFAULT NULL,
  `order` varchar(75) DEFAULT NULL,
  `family` varchar(75) DEFAULT NULL,
  `genus` varchar(75) DEFAULT NULL,
  `species` varchar(75) DEFAULT NULL,
  `organism` varchar(125) DEFAULT NULL,
  `sys_topHit` tinyint(1) NOT NULL DEFAULT '0',
  `fxn_topHit` tinyint(1) NOT NULL DEFAULT '0',
  `dateCreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `dateModified` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `deleted` tinyint(1) NOT NULL DEFAULT '0'
  );");
  $dbh->do( "CREATE INDEX blastn_query_name ON blastn(`query_name`);" );
  $dbh->do( "CREATE INDEX blastn_e_value ON blastn(`e_value`);" );
  $dbh->do( "CREATE INDEX blastn_sequenceId ON blastn(`sequenceId`);" );
  $dbh->do( "CREATE INDEX blastn_databasename ON blastn(`database_name`);" );
  $dbh->do( "CREATE INDEX blastn_sysTopHit ON blastn(`sys_topHit`);" );
  $dbh->do( "CREATE INDEX blastn_fxnTopHit ON blastn(`fxn_topHit`);" );
  v$dbh->do( "CREATE INDEX blastn_hitname ON blastn(`hit_name`);" );
  #$dbh->do( "CREATE INDEX blastn_query_name ON blastn(query_name,e_value);" );
  #$dbh->do( "CREATE INDEX blastn_domain ON blastn(`domain`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_kingdom ON blastn(`kingdom`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_phylum ON blastn(`phylum`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_class ON blastn(`class`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_order ON blastn(`order`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_family ON blastn(`family`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_genus ON blastn(`genus`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_species ON blastn(`species`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_organism ON blastn(`organism`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_hit_description ON blastn(`hit_description`,`deleted`,`e_value`);" );
  #$dbh->do( "CREATE INDEX blastn_qn_hn ON blastn(`query_name`,`hit_name`,`qry_start`,`hit_start`);" );


## Create table sequence
$dbh->do( "CREATE TABLE `sequence` (
  `id` INTEGER PRIMARY KEY AUTOINCREMENT,
  `libraryId` int(8) NOT NULL,
  `name` VARCHAR(255) NOT NULL,
  `header` TEXT NOT NULL,
  `gc` DOUBLE NOT NULL,
  `basepair` TEXT NOT NULL,
  `size` INT(8) NOT NULL,
  `comment` TEXT,
  `dateCreated` TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  `dateModified` TIMESTAMP,
  `deleted` tinyint(1) DEFAULT 0,
  `typeId` INT(8) NOT NULL
);");
$dbh->do( "CREATE INDEX sequence_libraryId ON sequence(`libraryId`);" );
#$dbh->do( "CREATE INDEX sequence_name ON sequence(`name`);" );

## Create table sequence_relationship
$dbh->do( "CREATE TABLE `sequence_relationship` (
  `subjectId` bigint(19) NOT NULL,
  `objectId` bigint(19) NOT NULL,
  `typeId` int(8) NOT NULL
);");
$dbh->do( "CREATE INDEX sequence_relationship_subjectId ON sequence_relationship(`subjectId`);" );
$dbh->do( "CREATE INDEX sequence_relationship_objectId ON sequence_relationship(`objectId`);" );
$dbh->do( "CREATE INDEX sequence_relationship_typeId ON sequence_relationship(`typeId`);" );

## Create table statistics
$dbh->do( "CREATE TABLE `statistics` (
  `id` INTEGER PRIMARY KEY AUTOINCREMENT,
  `libraryId` int(8) NOT NULL DEFAULT '0',
  `raw_read_cnt` int(11) DEFAULT NULL,
  `post_qc_read_cnt` int(11) DEFAULT NULL,
  `post_cdhit454_read_cnt` int(11) DEFAULT NULL,
  `read_cnt` int(11) NOT NULL DEFAULT '0',
  `read_mb` double NOT NULL DEFAULT '0',
  `complete_cnt` int(11) NOT NULL DEFAULT '0',
  `complete_mb` double NOT NULL DEFAULT '0',
  `complete_id` text,
  `incomplete_cnt` int(11) NOT NULL DEFAULT '0',
  `incomplete_mb` double NOT NULL DEFAULT '0',
  `incomplete_id` text,
  `lackstart_cnt` int(11) NOT NULL DEFAULT '0',
  `lackstart_mb` double NOT NULL DEFAULT '0',
  `lackstart_id` text,
  `lackstop_cnt` int(11) NOT NULL DEFAULT '0',
  `lackstop_mb` double NOT NULL DEFAULT '0',
  `lackstop_id` text,
  `archaea_cnt` int(11) NOT NULL DEFAULT '0',
  `archaea_mb` double NOT NULL DEFAULT '0',
  `archaea_id` text,
  `bacteria_cnt` int(11) NOT NULL DEFAULT '0',
  `bacteria_mb` double NOT NULL DEFAULT '0',
  `bacteria_id` text,
  `phage_cnt` int(11) NOT NULL DEFAULT '0',
  `phage_mb` double NOT NULL DEFAULT '0',
  `phage_id` text,
  `allviral_cnt` int(8) NOT NULL,
  `allviral_id` text NOT NULL,
  `topviral_cnt` int(8) NOT NULL,
  `topviral_id` text NOT NULL,
  `allmicrobial_cnt` int(8) NOT NULL,
  `allmicrobial_id` text NOT NULL,
  `topmicrobial_cnt` int(8) NOT NULL,
  `topmicrobial_id` text NOT NULL,
  `fxn_cnt` int(8) NOT NULL,
  `fxn_id` text,
  `unassignfxn_cnt` int(8) NOT NULL,
  `unassignfxn_id` text,
  `tRNA_cnt` int(8) NOT NULL DEFAULT '0',
  `tRNA_id` text,
  `rRNA_cnt` int(8) NOT NULL DEFAULT '0',
  `rRNA_id` text,
  `orfan_cnt` int(8) NOT NULL DEFAULT '0',
  `orfan_id` text,
  `lineage` text,
  `deleted` tinyint(4) NOT NULL DEFAULT '0',
  `dateCreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP
);");
#$dbh->do( "CREATE INDEX statistics_libraryId ON statistics(`libraryId`,`deleted`);" );
$dbh->do( "CREATE INDEX statistics_libraryId ON statistics(`libraryId`);" );

## Create table tRNA
$dbh->do( "CREATE TABLE `tRNA` (
  `id` INTEGER PRIMARY KEY AUTOINCREMENT,
  `sequenceId` int(8) NOT NULL,
  `num` tinyint(3) NOT NULL,
  `tRNA_start` mediumint(8) NOT NULL,
  `tRNA_end` mediumint(8) NOT NULL,
  `anti` varchar(6) NOT NULL,
  `intron` varchar(3) NOT NULL,
  `cove_start` mediumint(8) NOT NULL,
  `cove_end` mediumint(8) NOT NULL,
  `score` float NOT NULL,
  `dateCreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `deleted` tinyint(1) NOT NULL DEFAULT '0'
);");
$dbh->do( "CREATE INDEX tRNA_sequenceId ON tRNA(`sequenceId`);" );

## Create table mgol_library
$dbh->do( "CREATE TABLE `mgol_library` (
	`mgol_id` int(6) PRIMARY KEY,
	`ncbi_mg_src` varchar(100),
	`ncbi_parent_proj` varchar(50),
	`ncbi_proj` varchar(50),
	`ncbi_acc` varchar(100),
	`lib_name` varchar(50),
	`lib_prefix` varchar(100),
	`lib_shortname` varchar(50),
	`data_src` varchar(50),
	`url` varchar(200),
	`citation` text,
	`citation_pdf` varchar(100),
	`fasta_nt` varchar(100),
	`fasta_pep` varchar(100),
	`seq_type` varchar(100),
	`amplification` varchar(100),
	`seq_center` varchar(100),
	`seq_release_date` varchar(50),
	`lib_type` varchar(50),
	`na_type` varchar(50),
	`avg_read_len` float,
	`GC_pct` float,
	`filter_lower_um` varchar(50),
	`filter_upper_um` varchar(50),
	`virome` int(1),
	`project` varchar(100),
	`sample_date` varchar(20),
	`num_sites` int(4),
	`site_id` varchar(100),
	`genesis` varchar(100),
	`sphere` varchar(100),
	`ecosystem` varchar(100),
	`phys_subst` varchar(100),
	`extreme` int(1),
	`physio_chem_mods` varchar(100),
	`kg_climate` varchar(50),
	`host_tax` text,
	`host_common_name` varchar(100),
	`org_substr` varchar(200),
	`lat_zone` varchar(200),
	`depth_zone` varchar(200),
	`alt_zone` varchar(200),
	`region` varchar(200),
	`geog_place_name` varchar(200),
	`country` varchar(50),
	`lat_deg` float,
	`lat_hem` char(1),
	`lon_deg` float,
	`lon_hem` char(1),
	`Biome_EnvO` varchar(200),
	`Biome_desc` varchar(200),
	`Biome_EnvONum` varchar(100),
	`envo_feat_geo` varchar(100),
	`envo_feat_hab` varchar(100),
	`envo_feat_meso` varchar(100),
	`envo_matter` varchar(100),
	`envo_food` varchar(100),
	`env_annotator` varchar(100),
	`chl_ugkg` varchar(20),
	`chl_mgM3` varchar(20),
	`pres_atm` varchar(20),
	`sample_depth` varchar(20),
	`water_depth` varchar(20),
	`alt_m` varchar(20),
	`biomass_ugkg` varchar(20),
	`DIC_umolkg` varchar(20),
	`DIP_nmolkg` varchar(20),
	`DOC_umolkg` varchar(20),
	`DO_nmolkg` varchar(20),
	`NO3_nmolkg` varchar(20),
	`salinity_psu` varchar(20),
	`temp_c` varchar(20),
	`pH` varchar(20),
	`intel_prop` text,
	`qryDb` tinyint(1),
	`comments` text,
	`deleted` tinyint(1)
);");
$dbh->do( "CREATE INDEX lib_prefix ON mgol_library(`lib_prefix`);" );

open(OUT, ">", $options{pstore}."/mgol.sql") or $logger->logdie("Could not open file to write $options{pstore}/mgol.sql");
print OUT ".separator \"\\t\"\n";
print OUT ".import /opt/ergatis/autopipe_package/mgol_table.tab mgol_library";

close(OUT);

my $cmd = "sqlite3 $options{pstore}/processing.sqlite3";
$cmd .= " < $options{pstore}/mgol.sql";

system($cmd);

$logger->info("init database for complete");
exit(0);

###############################################################################
sub check_parameters {
    #### make sure sample_file and output_dir were passed
    my @required = qw(outdir pstore);

	foreach my $key (@required) {
		unless ($options{$key}) {
			pod2usage({-exitval => 2,  -message => "ERROR: Input $key not defined", -verbose => 1, -output => \*STDERR});
		    $logger->logdie("No input defined, plesae read perldoc $0\n\n");
        }
	}
}
