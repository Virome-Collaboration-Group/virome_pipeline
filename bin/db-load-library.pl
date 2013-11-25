#!/usr/bin/perl -w

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell


=head1 NAME

db-load-library.pl - load Library output to db

=head1 SYNOPSIS

USAGE: db-load-library.pl
            --id=/id the id of the library
            --user=/user/who/is/uploading/component
	    --outdir=/path/to/output/dir
	    --env=/env where to execute script
          [ --log=/path/to/logfile
            --debug=N
          ]

=head1 OPTIONS

B<--id, -i>
    The id number of the library you are running.
    You can find this in the library table on Jabba.
    
B<--user, -u>
    User who is uplaoding library infomation.

B<--outdir, o>
    Path to output dir where a tab file of library id and library prefix is stored.
    
B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script is used to load library infomation into MySQL database.

=head1  INPUT

The input to this is defined using the --id.  This should be
the id number of the library you wish to run. If the library you
wish to run does not have an id than you will need to add
the information tot he table.

=head1  CONTACT

    Daniel J. Nasko
    dan.nasko@gmail.com

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
BEGIN {
  use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
                          'id|i=s',
                          'user|u=s',
			  'outdir|o=s',
			  'loc|e=s',
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

#############################################################################
#### DEFINE GLOBAL VAIRABLES.
##############################################################################
my $db_user;
my $db_pass;
my $dbname;
my $db_host;
my $host;
my @RESULTS;
my $user = $options{user};
my $id = $options{id};


my $dbh;

print STDOUT ("CHECKING PARAMS\n");

## make sure everything passed was peachy
&check_parameters(\%options);
##############################################################################

print STDOUT ("Checking production server for library information\n");

## MySQL statement to gather library information
my $lib_info_sql = qq/SELECT id, name, prefix, server FROM library WHERE id = ? and user like "\%"  ? "\%" and progress = "standby"/;
my $edit_lib_sql = qq/UPDATE library SET progress = "running" where id = ?/;
  
## prep statement
my $sth_info = $dbh->prepare($lib_info_sql);
my $sth_edit = $dbh->prepare($edit_lib_sql);

print STDOUT ("PROCESS FILE\n");

#execute to check lib existance.
$sth_info->execute($id, $user);
$sth_edit->execute($id);

#raise error.
$sth_info->{RaiseError} = 1;

#$outfile =~ s/\s+/_/g;
my $outfile = $options{outdir}."/"."db-load-library.txt";
    
open (LIB, ">$outfile") or logger->logdie("Cannot open output file $outfile\n");

while (@RESULTS = $sth_info->fetchrow_array) {	print LIB "$RESULTS[0]\t$RESULTS[1]\t$RESULTS[2]\t$RESULTS[3]";}
print LIB "\t", $options{'loc'};
close LIB;

$dbh->disconnect();
close(DAT);
exit(0);

##############################################################################
# SUBS
##############################################################################
sub check_parameters {
    
    ## at least one input type is required
    unless ( $options{id} && $options{user} && $options{outdir} && $options{loc}) {
	$logger->logdie("No input defined, plesae read perldoc $0\n");
        exit(1);
    }

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
    }
    
    if ($options{loc} eq 'dbi'){
	$db_user = q|bhavsar|;
	$db_pass = q|P3^seus|;
	$dbname = q|VIROME|;
	$db_host = q|virome.dbi.udel.edu|;
	$host = q|virome.dbi.udel.edu|;
    }elsif ($options{loc} eq 'diag'){
        $db_user = q|dnasko|;
        $db_pass = q|dnas_76|;
        $dbname = q|virome|;
        $db_host = q|virome-db.igs.umaryland.edu|;
        $host = q|virome-db.igs.umaryland.edu|;
    }elsif ($options{loc} eq 'diag1'){
        $db_user = q|dnasko|;
        $db_pass = q|dnas_76|;
        $dbname = q|virome|;
        $db_host = q|virome-db.igs.umaryland.edu|;
        $host = q|virome-db.igs.umaryland.edu|;
    }elsif ($options{loc} eq 'diag2'){
        $db_user = q|dnasko|;
        $db_pass = q|dnas_76|;
        $dbname = q|virome|;
        $db_host = q|virome-db.igs.umaryland.edu|;
        $host = q|virome-db.igs.umaryland.edu|;
    }elsif ($options{loc} eq 'diag3'){
        $db_user = q|dnasko|;
        $db_pass = q|dnas_76|;
        $dbname = q|virome|;
        $db_host = q|virome-db.igs.umaryland.edu|;
        $host = q|virome-db.igs.umaryland.edu|;
    }elsif ($options{loc} eq 'diag4'){
        $db_user = q|dnasko|;
        $db_pass = q|dnas_76|;
        $dbname = q|virome|;
        $db_host = q|virome-db.igs.umaryland.edu|;
        $host = q|virome-db.igs.umaryland.edu|;
    }elsif ($options{loc} eq 'diag5'){
        $db_user = q|dnasko|;
        $db_pass = q|dnas_76|;
        $dbname = q|virome|;
        $db_host = q|virome-db.igs.umaryland.edu|;
        $host = q|virome-db.igs.umaryland.edu|;
    }elsif ($options{loc} eq 'igs'){
	$db_user = q|dnasko|;
	$db_pass = q|dnas_76|;
	$dbname = q|virome|;
	$db_host = q|jabba.igs.umaryland.edu|;
	$host = q|jabba.igs.umaryland.edu|;
    }elsif ($options{loc} eq 'camera'){
        $db_user = q|virome_app|;
        $db_pass = q|camera123|;
        $dbname = q|virome_stage|;
        $db_host = q|dory.internal.crbs.ucsd.edu|;
        $host = q|dory.internal.crbs.ucsd.edu|;
    }elsif ($options{loc} eq 'ageek') {
	$db_user = q|bhavsar|;
	$db_pass = q|Application99|;
	$dbname = q|virome|;
	$db_host = q|10.254.0.1|;
	$host = q|10.254.0.1|;
    }else {
	$db_user = q|kingquattro|;
	$db_pass = q|Un!c0rn|;
	$dbname = q|VIROME|;
	$db_host = q|localhost|;
	$host = q|localhost|;
    }
    
    $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$db_host",
	"$db_user", "$db_pass",{PrintError=>1, RaiseError =>1, AutoCommit =>1}) || die print STDERR("Cound not open db connection\n");
}

##############################################################################
sub getServer{
  my $env = $_[0];
  my @server = ("calliope", "polyhymnia", "terpsichore", "thalia");
  
  if ($env =~ /water/i){
    return $server[2];
  } elsif ($env =~ /soil/i){
    return $server[1];
  } elsif ($env =~ /extreme/i){
    return $server[3];
  } elsif ($env =~ /solid substrate/i){
    return $server[3];
  } else {
    return $server[0];
  }
}

##############################################################################
sub trim($) {
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}
