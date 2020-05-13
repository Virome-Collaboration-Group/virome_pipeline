#!/usr/bin/perl

=head1 NAME

virome_little_run_pipeline.pl - A script allowing VIROME-Ergatis pipeline instantiation
and execution via the API. It's only going to do so for a little virome pipeline though.
i.e. 2 components.

=head1 SYNOPSIS

USAGE: virome_little_run_pipeline.pl
            --template_directory=/path/to/some_dir/
            --repository_root=/path/to/project_dir
            --id_repository=/path/to/foo/id_repository
            --ergatis_ini=/path/to/ergatis.ini
            --fasta=/path/to/file.fasta
            --threads=number_of_max_threads_to_use

=head1 OPTIONS

B<--template_directory,-t>
    The full path to a directory containing a pipeline.layout file and its associated
    INI configs

B<--repository_root,-r>
    Also known as a project area, under this directory we should find

B<--id_repository,-i>
    Path to some global-level ID repository.  Will only be used to pull a pipeline ID.

B<--ergatis_ini,-e>
    Path to the ergatis.ini file, usually kept in the cgi directory of an interface instance.

B<--fasta,-f>
    Path to the input FASTA file

B<--threads,-d>
    Number of max threads to use during execution of workflow

B<--log,-l>
    Log file

B<--help,-h>
    This help message


=head1  DESCRIPTION

It's not uncommon to want to instantiate and launch a pipeline from the command-line rather
than using the web interface.  This script illustrates how to do that in just a few lines of
code (4 lines, really, the rest is common perl script template code).

=head1  INPUT

Described in the options above, but the main input is the directory containing an Ergatis
pipeline template.  Once the template is loaded, you can use the Ergatis::ConfigFiles module
to customize the component configuration files for your application.  An example is shown in the
code.

=head1  OUTPUT

This script will use the pipeline template to instantiate a full pipeline in the passed project
area (repository root) and then execute it.

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

    Daniel Nasko
    dan.nasko@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use FindBin;
use File::Path qw(make_path remove_tree mkpath);
use JSON qw(encode_json decode_json);
use lib "$FindBin::Bin/../lib";
use lib ("/var/www/html/cgi") ;

use Ergatis::ConfigFile;
use Ergatis::SavedPipeline;

END { _log("The program ran for " . (time() - $^T) . " seconds"); }

my %options = ();
my $results = GetOptions (\%options,
                          'fasta|f=s',
                          'template_directory|t=s',
                          'repository_root|r=s',
                          'id_repository|i=s',
                          'ergatis_ini|e=s',
                          'threads|d=s',
                          'log|l=s',
                          'debug|v=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

print "DEBUG VALUE: " .$options{debug}. "\n";

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
# if (defined $options{log}) {
#     open($logfh, ">$options{log}") || die "can't create log file: $!";
# }

##############################
## important bits here
##############################
my $fasta = $options{'fasta'};
my $cmd = "";

my $filename = fileparse($fasta, qr/\.[^.]*/);
my $output_dir = "/opt/output/". $filename ."_". timestamp();

print "DEBUG VALUE: " .$options{debug}. "\n";

&create_output_dir($output_dir);

my $version_info = parse_version_info();

my $template = Ergatis::SavedPipeline->new(
               template => "$options{template_directory}/pipeline.blastonly.layout");

my $pipeline = $template->write_pipeline( repository_root => $options{repository_root},
                                          id_repository => $options{id_repository} );

open($logfh, ">", "$output_dir/logs/$pipeline->{'id'}.analytic") or
    die "Could not open file to write error log $output_dir/logs/$pipeline->{'id'}.analytic\n";

_log("VIROME Analysis Pipeline analytics.");
_log("Start time: " . format_time());
_log("Version info: " . $version_info->{version} . " build date: " . $version_info->{date});
_log("Input file name: " . $fasta);
_log("Input file size: " . sprintf("%.2f", ((-s $fasta)/(1024 * 1024))) . "MB");

my $no_of_seq = `grep -c "^>" $fasta`;
chomp $no_of_seq;

_log("Number of raw input seq: " . $no_of_seq);

## here you can use Ergatis::ConfigFiles to edit some of the newly-written
##  component configurations before pipeline execution.  One example is shown.
##  naming and path conventions allow you to know where the component file is

#### set output dir for init-db
my $init_db_config = new Ergatis::ConfigFile(
    -file => "$options{repository_root}/workflow/runtime/init-db/" . $pipeline->id . "_default/init-db.default.user.config");
$init_db_config->setval('parameters', '$;PERSISTENT_STORAGE$;', $output_dir );
$init_db_config->RewriteConfig();

#### final dump input file setup
my $dump_db_config = new Ergatis::ConfigFile(
    -file => "$options{repository_root}/workflow/runtime/results_blastonly/" . $pipeline->id . "_default/results_blastonly.default.user.config");
$dump_db_config->setval('input', '$;INPUT_FILE$;', $fasta);
$dump_db_config->setval('parameters', '$;PERSISTENT_STORAGE$;', $output_dir);
$dump_db_config->setval('parameters', '$;VERBOSE63$;', $options{debug});
$dump_db_config->setval('parameters', '$;DATABASE_FILE$;', $output_dir . '/processing.sqlite3' );
$dump_db_config->RewriteConfig();

#### $fasta_size_filter
my $fasta_size_filter_config = new Ergatis::ConfigFile(
    -file => "$options{repository_root}/workflow/runtime/fasta_size_filter/" . $pipeline->id . "_default/fasta_size_filter.default.user.config");
$fasta_size_filter_config->setval('input', '$;INPUT_FILE$;', $fasta );
$fasta_size_filter_config->RewriteConfig();

#### set univec subject database name and path
my $univec_config = new Ergatis::ConfigFile(
    -file => "$options{repository_root}/workflow/runtime/ncbi-blastn-plus/" . $pipeline->id . "_univec/ncbi-blastn-plus.univec.user.config");
$univec_config->setval('parameters', '$;DATABASE_PATH$;', '/opt/database/' . $version_info->{univec});
$univec_config->RewriteConfig();

#### set rRNA subject database name and path
my $rna_config = new Ergatis::ConfigFile(
    -file => "$options{repository_root}/workflow/runtime/ncbi-blastn-plus/" . $pipeline->id . "_rna/ncbi-blastn-plus.rna.user.config");
$rna_config->setval('parameters', '$;DATABASE_PATH$;', '/opt/database/' . $version_info->{rna});
$rna_config->RewriteConfig();

#### set max threads limit for rubble blast, subject database name and database path
my $uniref_rubble_config = new Ergatis::ConfigFile(
    -file => "$options{repository_root}/workflow/runtime/rubble/" . $pipeline->id . "_uniref/rubble.uniref.user.config");
$uniref_rubble_config->setval('parameters', '$;THREADS$;', $options{threads} );
$uniref_rubble_config->setval('parameters', '$;DATABASE_PATH$;', '/opt/database/' . $version_info->{uniref} );
$uniref_rubble_config->setval('parameters', '$;DATABASE_CLUST_PATH$;', '/opt/database/' . $version_info->{uniref_clust} );
$uniref_rubble_config->setval('parameters', '$;LOOKUP$;', '/opt/database/' . $version_info->{uniref_lookup} );
$uniref_rubble_config->RewriteConfig();

my $mgol_rubble_config = new Ergatis::ConfigFile(
    -file => "$options{repository_root}/workflow/runtime/rubble/" . $pipeline->id . "_mgol/rubble.mgol.user.config");
$mgol_rubble_config->setval('parameters', '$;THREADS$;', $options{threads} );
$mgol_rubble_config->setval('parameters', '$;DATABASE_PATH$;', '/opt/database/' . $version_info->{mgol} );
$mgol_rubble_config->setval('parameters', '$;DATABASE_CLUST_PATH$;', '/opt/database/' . $version_info->{mgol_clust} );
$mgol_rubble_config->setval('parameters', '$;LOOKUP$;', '/opt/database/' . $version_info->{mgol_lookup} );
$mgol_rubble_config->RewriteConfig();

my $tRNAscan_config = new Ergatis::ConfigFile(
    -file => "$options{repository_root}/workflow/runtime/tRNAscan-SE/" . $pipeline->id . "_default/tRNAscan-SE.default.user.config");
$tRNAscan_config->setval('parameters', '$;THREADS$;', $options{threads} );
$tRNAscan_config->RewriteConfig();


#### point to PERSISTENT_STORAGE sqlite3 file
my @array = qw(mgol rna uniref univec);
foreach my $token (@array) {
    my $pstore_config = new Ergatis::ConfigFile(
        -file => "$options{repository_root}/workflow/runtime/blast-result-prep/" . $pipeline->id . "_${token}/blast-result-prep.${token}.user.config");
    $pstore_config->setval('parameters', '$;DATABASE_FILE$;', $output_dir . '/processing.sqlite3' );
    $pstore_config->RewriteConfig();
}

#### clean_expand_btab sqlite3 file location
foreach my $token (@array) {
    my $pstore_config = new Ergatis::ConfigFile(
        -file => "$options{repository_root}/workflow/runtime/clean_expand_btab/" . $pipeline->id . "_${token}/clean_expand_btab.${token}.user.config");
    $pstore_config->setval('parameters', '$;DATABASE_FILE$;', $output_dir . '/processing.sqlite3' );
    $pstore_config->RewriteConfig();
}

#### stats script sqlite3 file location
# undef @array;
# @array = qw(env_lib_stats fxnal_bin_all_db fxnal_bin_per_db gen_lib_stats libraryHistogram orfan sequence_classification taxonomy_binning);
# foreach my $component (@array) {
#     my $stats_config = new Ergatis::ConfigFile(
#         -file => "$options{repository_root}/workflow/runtime/${component}/" . $pipeline->id . "_default/${component}.default.user.config");
#     $stats_config->setval('input', '$;INPUT_FILE$;', $output_dir . '/processing.sqlite3' );
#     $stats_config->RewriteConfig();
# }

#### db_upload script sqlite3 file location
undef @array;
@array = qw(orf_nuc orf_pep rna-blast sequence_read sequence_relationship sequence_rna trna univec-blast);
foreach my $token (@array) {
    my $pstore_config = new Ergatis::ConfigFile(
        -file => "$options{repository_root}/workflow/runtime/db-upload/" . $pipeline->id . "_${token}/db-upload.${token}.user.config");
    $pstore_config->setval('parameters', '$;DATABASE_FILE$;', $output_dir . '/processing.sqlite3' );
    $pstore_config->RewriteConfig();
}

#### sequence-prep sqlite3 file location
undef @array;
@array = qw(orf_pep orf_nuc rna-clean rna);
foreach my $token (@array) {
    my $pstore_config = new Ergatis::ConfigFile(
        -file => "$options{repository_root}/workflow/runtime/sequence-prep/" . $pipeline->id . "_${token}/sequence-prep.${token}.user.config");
    $pstore_config->setval('parameters', '$;DATABASE_FILE$;', $output_dir . '/processing.sqlite3' );
    $pstore_config->RewriteConfig();
}

#### database dump sqlite3 location
undef @array;
@array = qw(tRNAscan-prep);
foreach my $component (@array) {
    my $pstore_config = new Ergatis::ConfigFile(
        -file => "$options{repository_root}/workflow/runtime/${component}/" . $pipeline->id . "_default/${component}.default.user.config");
    $pstore_config->setval('parameters', '$;DATABASE_FILE$;', $output_dir . '/processing.sqlite3' );
    $pstore_config->RewriteConfig();
}

#### sequence_relationship-prep sqlite3 location
my $pstore_config = new Ergatis::ConfigFile(
    -file => "$options{repository_root}/workflow/runtime/sequence_relationship-prep/" . $pipeline->id . "_default/sequence_relationship-prep.default.user.config");
$pstore_config->setval('input', '$;INPUT_FILE$;', $output_dir . '/processing.sqlite3' );
$pstore_config->RewriteConfig();


## Get ready to rumble . . .
my $ergatis_cfg = new Ergatis::ConfigFile( -file => $options{ergatis_ini} );

# The 'block' term is used to make the pipeline invocation synchronous, so that
# we can determine success/failure by the return and use that to determine this
# script's exit value.
my $success = $pipeline->run(ergatis_cfg => $ergatis_cfg, block => 1);

if (! $success) {
    my $stderr = "Problem running pipeline id:$pipeline->{'id'} with input: $options{'fasta'}\n\n";
    $stderr .= "$pipeline->{'diagnostics'}->{'complete_components'} of ";
    $stderr .= "$pipeline->{'diagnostics'}->{'total_components'} completed\n";
    $stderr .= "\n";

    $stderr .= "ERROR running component(s): \n";
    foreach my $c (@{$pipeline->{'diagnostics'}->{'components'}}) {
        my $t = $c;
        $t =~ s/\.default$//;

        $stderr .= "\t$c\n";

        foreach my $t (@{$pipeline->{'diagnostics'}->{'command_info'}->{$c}}) {
            my $c = $$t[0];
            my $f = $$t[1];

            if (length $f) {
                open(FHD, "<", $f) or die "Could not open file $f\n$!";
                while(<FHD>) {
                    $stderr .= "\t\t$_";
                }
                $stderr .= "\n";
            }
        }
    }

    print STDERR $stderr;

    $cmd = "scp $options{repository_root}/workflow/runtime/pipeline/$pipeline->{'id'}/pipeline.xml.log";
    $cmd .= " $output_dir/logs/.";
    system("$cmd");

    open(OUT, ">", "$output_dir/logs/$pipeline->{'id'}.stderr") or
        die "Could not open file to write error log $output_dir/logs/$pipeline->{'id'}.stderr\n";
    print OUT $stderr;
    close(OUT);
} else {
    #### pipeline ended in success collect few more data points.

    #### extract number of orfs
    my $dir = "$options{repository_root}/output_repository/mga2seq_pep/1_default/i1/g1";

    opendir(DIR, "$dir") or die "Could open directry $dir";
    while( ($filename = readdir(DIR))) {
       if ($filename =~ /\.pep$/) {
           my $no_of_orf = `grep -c "^>" $dir/$filename`;
           chomp $no_of_orf;

           _log("Number of ORFs: " . $no_of_orf);
           my $v = "$dir/$filename";
           _log("File size of ORFs: " . sprintf("%.2f", ((-s $v)/(1024 * 1024))) . "MB");
       }
    }
    closedir(DIR);

    #### read and write read/contig stats
    $dir = "$options{repository_root}/output_repository/nt_fasta_check/1_default/i1/g1";

    opendir(DIR, "$dir") or die "Could open directry $dir";
    while( ($filename = readdir(DIR))) {
       if ($filename =~ /\.stats$/) {
           open(S, "<", "$dir/$filename") or die "Could not open file $dir/$filename\n";
           while(<S>) {
               _log($_);
           }
           close(S);
       }
    }
    closedir(DIR);

    #### blastp uniref output
    $dir = "$options{repository_root}/output_repository/clean_expand_btab/1_uniref/i1/g1";

    opendir(DIR, "$dir") or die "Could open directry $dir";
    while( ($filename = readdir(DIR))) {
       if ($filename =~ /\.btab/) {
           my $len = `wc -l $dir/$filename`;
           chomp $len;

           _log("Number of BLASTP rows (UNIREF100P): " . $len . "\n");
       }
    }
    closedir(DIR);

    $dir = "$options{repository_root}/output_repository/clean_expand_btab/1_mgol/i1/g1";

    opendir(DIR, "$dir") or die "Could open directry $dir";
    while( ($filename = readdir(DIR))) {
       if ($filename =~ /\.btab/) {
           my $len = `wc -l $dir/$filename`;
           chomp $len;

           _log("Number of BLASTP rows (MGOL): " . $len . "\n");
       }
    }
    closedir(DIR);
}

#print Dumper(\$pipeline->{'diagnostics'});

## Determine the exit value based on the success of the pipeline.
my $exit_value = ($success) ? 0 : 1;
_log("End time: " . format_time());

close(LOG);

exit $exit_value;

###############################################################################
sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

###############################################################################
sub check_parameters {
    #### make sure required arguments were passed
    my @required = qw( template_directory repository_root id_repository ergatis_ini fasta);

    foreach my $key (@required) {
        unless ($options{$key}) {
            die "--$key is a required option";
        }
    }

    $options{debug} = 0 unless (defined $options{debug});
}
###############################################################################
sub create_output_dir {
    my $dir = shift;

    #mkpath($dir.'/logs', 0, '0755');
    make_path($dir.'/logs');
}

###############################################################################
sub parse_version_info {
    my $data;
    {
        local $/ = undef;
        open (FHD, "<", '/opt/database/version.json.current')
            or die "Could not open file /opt/database/version.json.current to read\n";
        $data = <FHD>;
        close FHD;
    }

    return decode_json($data);
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

###############################################################################
sub format_time {

    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();

    my @months   = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

    my $year    = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";

    return $theTime;
}
###############################################################################
