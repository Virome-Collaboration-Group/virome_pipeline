#!/usr/bin/perl

=head1 NAME
   FxnalDBBreak.pl

=head1 SYNOPSIS

    USAGE: FxnalDBBreak.pl --server server-name --env dbi [--library libraryId]

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
    Create XML document that contaions information to draw Fxnal DB
    drill down chart.

=head1  INPUT
    The input is defined with --server,  --library.

=head1  OUTPUT
   Updated blastp tables for all/specifed library.

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   FxnalDBBreak.pl --server calliope --env dbi --library 31

=cut

use strict;
use warnings;
use IO::File;
use DBI;
use UTILS_V;
use Tree::Nary;
use Pod::Usage;
use Data::Dumper;
use JSON;
use File::Path qw(make_path remove_tree mkpath);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

BEGIN {
  use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
            			  'input|i=s',
            			  'outdir|o=s',
                          'lookupdb|b=s',
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
my $this;
my $cmd = "";

$this->{output_dir} = "$options{outdir}/xDocs";

my $libId = 1;

my $blst_stmt = qq{SELECT	b.hit_name, b.sequenceId, b.database_name, s.libraryId
				   FROM		blastp b INNER JOIN sequence s
				     ON 	s.id=b.sequenceId
				   WHERE	s.libraryId = ?
						AND b.e_value <= 0.001
						AND b.fxn_topHit = 1
						AND b.database_name = ?}; #change kegg to ?

my @dbArray = ("ACLAME", "COG", "UNIREF100P", "KEGG", "SEED", "PHGSEED");

my $r_node;
my $subject_db_hash = {};
$r_node->{name} = "root";
$r_node->{value} = 1;
$r_node->{id_list} = "NA";

my $tree = Tree::Nary->new($r_node);

#### process list for each functional database
foreach my $db (@dbArray){
	print "\tProcessing db $db\n";

    #### Load subjec data into memory so subsequent lookup are faster.
    load_subject_db($db);

	#### create new tree
	$tree = Tree::Nary->new($r_node);
	my $blst_qry = $dbh->prepare($blst_stmt);
	$blst_qry->execute($libId, $db);

	#### for each database process top fxn hit per sequence
	while (my $hits = $blst_qry->fetchrow_hashref()) {
		my $local_tree;
		if ($$hits{database_name} =~ /uniref100p|aclame/i){
			$local_tree = createLocalTreeAU($$hits{hit_name}, $$hits{sequenceId}, $$hits{database_name});
		} else {
			$local_tree = createLocalTreeR($$hits{hit_name}, $$hits{sequenceId}, $$hits{database_name});
		}
		#$tree = fuse_tree($tree,$local_tree);
        fuse_tree($local_tree);
	}

	my $xDoc;
	my $iDoc;
	my $temp;
	($xDoc, $iDoc, $temp) = node_to_xml($tree, 0, 1, $db, $libId);

	my $n = uc $db;
	if ($n =~ /uniref100p/i){
		$n = "GO";
	}

	print "\t\tWriting xDoc file\n";
	open (FHD, ">", $this->{output_dir} ."/". $n."_XMLDOC_".$libId.".xml") or $logger->logdie("Can not write to file ".$this->{output_dir} ."/". $n."_XMLDOC_".$libId.".xml");
	print FHD "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	print FHD "<root>\n";
	print FHD $xDoc;
	print FHD "</root>";
	close FHD;

	print "\t\tWriting idDoc file\n";
	open (FHD, ">", $this->{output_dir} ."/". $n."_IDDOC_".$libId.".xml") or $logger->logdie("Can not write to file ".$this->{output_dir} ."/". $n."_IDDOC_".$libId.".xml");
	print FHD "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	print FHD "<root>\n";
	print FHD $iDoc;
	print FHD "</root>";
	close FHD;
}

$dbh->disconnect;


$logger->info("Functional binning per database complete");
exit(0);


###############################################################################
####  SUBS
###############################################################################
sub check_parameters {
    my $options = shift;

    my @required = qw(input outdir lookupdb);

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

###############################################################################
sub load_subject_db {
    my $db = shift;

    my $dbh1 = DBI->connect("dbi:SQLite:dbname=$options{lookupdb}",
                            "", "", { RaiseError => 1}
                            ) or die $DBI::errstr;

    $db = lc $db;
    $db =~ s/uniref100p/go/;  #### change table name from uniref100p to go

    $subject_db_hash = {};
    my $acc_stmt = "";
    $acc_stmt = qq{SELECT realacc, block FROM ${db}_lookup};

    my $acc_qry = $dbh1->prepare($acc_stmt);
    $acc_qry->execute();
    $subject_db_hash = $acc_qry->fetchall_hashref("realacc");
}
###############################################################################
sub createLocalTreeR {
	my ($hit_name, $sequenceId, $database) = (shift,shift,shift);
	my $local_tree = Tree::Nary->new($r_node);

	my @acc = split(/;/,$hit_name);

	my $node = $local_tree->{children};

    #### get the array of hash for the give accession
    #### from subject_db_hash
    my $block = $subject_db_hash->{$acc[0]}->{block};

    #### end fuction if no annoation block for give acc
    #### e.g go acc A0A0F9WK84
    return $local_tree if (length($block) == 0);

    foreach my $anno ( @{ decode_json($block) } ) {
        my $n_node = node_child_exist($local_tree, $anno->{fxn1});
        if (! defined $n_node){
            my $new_node;
            $new_node->{name} = $anno->{fxn1};
            $new_node->{value} = 1;
            $new_node->{id_list} = $sequenceId;
            $node = Tree::Nary->append_data($local_tree, $new_node);
        } else {
            $node = $n_node;
        }

		$n_node = node_child_exist($node, $anno->{fxn2});
		if (! defined $n_node){
			my $new_node;
			$new_node->{name} = $anno->{fxn2};
			$new_node->{value} = 1;
			$new_node->{id_list} = $sequenceId;
			$node = Tree::Nary->append_data($node, $new_node);
		} else {
			$node = $n_node;
		}

		$n_node = node_child_exist($node, $anno->{fxn3});
		if (! defined $n_node){
			my $new_node;
			$new_node->{name} = $anno->{fxn3};
			$new_node->{value} = 1;
			$new_node->{id_list} = $sequenceId;
			$node = Tree::Nary->append_data($node, $new_node);
		} else {
			$node = $n_node;
		}
	}

	return $local_tree;
}

###############################################################################
sub createLocalTreeAU {
	my ($hit_name, $sequenceId, $database) = (shift,shift,shift);
	my $local_tree = Tree::Nary->new($r_node);
	my $acc_stmt = "";

	my @acc = split(/;/, $hit_name);

	my $node = $local_tree->{children};
	my $prev_chain_id = -1;

    #### get the array of hash for the give accession
    #### from subject_db_hash
    my $block = $subject_db_hash->{$acc[0]}->{block};

    #### end fuction if no annoation block for give acc
    #### e.g go acc A0A0F9WK84
    return $local_tree if (length($block) == 0);

    foreach my $anno ( @{ decode_json($block) } ) {
    	my $n_node;
		if ($prev_chain_id != $anno->{chain_id}){
			$n_node = node_child_exist($local_tree, $anno->{name});
		} else {
			$n_node = node_child_exist($node, $anno->{name});
		}

		if (!defined $n_node){
			my $new_node;
			$new_node->{name} = $anno->{name};
			$new_node->{value} = 1;
			$new_node->{id_list} = $sequenceId;

			if ($prev_chain_id != $anno->{chain_id}){
				$node = Tree::Nary->append_data($local_tree, $new_node);
			} else {
				$node = Tree::Nary->append_data($node, $new_node);
			}
		} else {
			$node = $n_node;
		}

		$prev_chain_id = $anno->{chain_id};
	}

	return $local_tree;
}

###############################################################################
sub node_child_exist {
	my ($node, $name) = (shift, shift);

	my $t_node = $node->{children};

	while (defined $t_node) {
		if (lc($t_node->{data}{name}) eq lc($name)){
			return $t_node;
		}
		$t_node = $t_node->{next};
	}

	return undef;
}

###############################################################################
sub fuse_tree {
	#my ($tree, $local) = (shift, shift);
    my ($local) = (shift);

	if (Tree::Nary->is_root($tree) && Tree::Nary->is_leaf($tree)){
		return $local;
	}

	if ((!defined $tree) && (defined $local)){
		return $local;
	}

	if (!defined $local || !defined $tree){
		return $tree;
	}

	for (my $i=0; $i < Tree::Nary->n_children($local); $i++){
		my $flag = 1;
		for (my $j=0; $j < Tree::Nary->n_children($tree); $j++){
			if (Tree::Nary->nth_child($tree,$j)->{data}{name} eq Tree::Nary->nth_child($local,$i)->{data}{name}){
				Tree::Nary->nth_child($tree,$j)->{data}{value}++;
				Tree::Nary->nth_child($tree,$j)->{data}{id_list}.=",".Tree::Nary->nth_child($local,$i)->{data}{id_list};
				$flag = 0;
				fuse_tree(Tree::Nary->nth_child($tree,$j),Tree::Nary->nth_child($local,$i));
			}
		}

		if ($flag){
			my $node = Tree::Nary->append_data($tree,Tree::Nary->nth_child($local,$i)->{data});
			fuse_tree($node,Tree::Nary->nth_child($local,$i));
		}
	}

	#return $tree;
}

###############################################################################
sub node_list {
	my ($node, $ref_arg) = 	(shift, shift);

	if (defined($ref_arg)){
		$$ref_arg .= $node->{data}{name}."[".$node->{data}{value}."] [".$node->{data}{id_list}."]\n";
	}

	return($Tree::Nary::FALSE);
}

###############################################################################
sub node_to_xml {
	my ($node,$level,$tag,$database,$id) = (shift,shift,shift,shift,shift);

	my $db = ($database =~ /UNIREF100P/) ? "GO" : $database;

	if (!defined $node){
		return ("","",$tag);
	}

	my $tab = "";
	my $xDoc = "";
	my $iDoc = "";

	for (my $i=0; $i<$level; $i++){
		$tab .= "\t";
	}

	for (my $i=0; $i < Tree::Nary->n_children($node); $i++){
		$xDoc .= $tab. "<FUNCTION_" .(Tree::Nary->depth(Tree::Nary->nth_child($node,$i))-1);
		$xDoc .= " NAME=\"" .Tree::Nary->nth_child($node,$i)->{data}{name}. "\" LABEL=\"";
		$xDoc .= Tree::Nary->nth_child($node,$i)->{data}{name}. "\" VALUE=\"";
		$xDoc .= Tree::Nary->nth_child($node,$i)->{data}{value}. "\" TAG=\"TAG_" .$tag. "\"";
		$xDoc .= " IDFNAME=\"" .$db."_IDDOC_$id.xml\"";

		$iDoc .= "<TAG_" .$tag++. " IDLIST=\"" .Tree::Nary->nth_child($node,$i)->{data}{id_list}. "\"/>\n";

		(my $xd, my $id, $tag) = node_to_xml(Tree::Nary->nth_child($node,$i),$level+1,$tag,$database,$id);
		if (length($xd) > 0){
			$xDoc .= ">\n"; #close pevious tag;
			$xDoc .= $xd;
			$xDoc .= $tab . "</FUNCTION_" .(Tree::Nary->depth(Tree::Nary->nth_child($node,$i))-1). ">\n";
			$iDoc .= $id;
		} else {
			$xDoc .= "/>\n";
		}
	}

	return ($xDoc, $iDoc, $tag);
}

###############################################################################
sub make_proper_case {
	my ($string) = @_;
	my @words = split (/\s+/, lc $string);
	my @new_words = ();
	my $new_word = "";

	foreach my $word (@words) {
		# Starts with Non-Alphanum Character
		my $starting_non_alphanum = "";
		if ($word =~ /^(\W)+(.*)/) {
			$starting_non_alphanum = $1;
			$word = $2;
		}

		# Ends with Non-Alphanum Character
		my $ending_non_alphanum = "";
		if ($word =~ /(.*)(\W)+$/) {
			$word = $1;
			$ending_non_alphanum = $2;
		}

		# Contains a Non-Alphanum Character
		if ($word =~ /^(\w+)(\W)(\w+)(\W?)(\w?)$/) {
			my $p1_word = $1;
			my $p2_non_alphanum = $2;
			my $p3_word = $3;
			my $p4_non_alphanum = $4;
			my $p5_letter = $5;

			$p1_word = ucfirst $p1_word;
			$p5_letter = lc $p5_letter;

			if (length $p1_word > 2 && length $p3_word == 1) {
				$p3_word = lc $p3_word;
			} elsif (length $p1_word == 1 && length $p3_word == 1) {
				$p3_word = uc $p3_word;
			} else {
				$p3_word = ucfirst $p3_word;
			}

			$new_word = $p1_word . $p2_non_alphanum . $p3_word . $p4_non_alphanum . $p5_letter;
		# Other
		} else {
			$new_word = ucfirst $word;
		}

		# Recombine the Alphanum Character
		$new_word = $starting_non_alphanum . $new_word . $ending_non_alphanum;

		push (@new_words, $new_word);
	}
	my $new_string = join(" ", @new_words);

	$new_string =~ s/(\w,?) And (\w)/$1 and $2/g;
	$new_string =~ s/(\w,?) Or (\w)/$1 or $2/g;
	$new_string =~ s/(\w,?) But (\w)/$1 but $2/g;

	$new_string =~ s/(\w) At (\w)/$1 at $2/g;
	$new_string =~ s/(\w) In (\w)/$1 in $2/g;
	$new_string =~ s/(\w) On (\w)/$1 on $2/g;
	$new_string =~ s/(\w) To (\w)/$1 to $2/g;
	$new_string =~ s/(\w) From (\w)/$1 from $2/g;

	$new_string =~ s/(\w) Is (\w)/$1 is $2/g;
	$new_string =~ s/(\w) A (\w)/$1 a $2/g;
	$new_string =~ s/(\w) An (\w)/$1 an $2/g;
	$new_string =~ s/(\w) Am (\w)/$1 am $2/g;
	$new_string =~ s/(\w) For (\w)/$1 for $2/g;
	$new_string =~ s/(\w) Of (\w)/$1 of $2/g;
	$new_string =~ s/(\w) The (\w)/$1 the $2/g;

	if (length $new_string > 60) {
		$new_string =~ s/(\w) With (\w)/$1 with $2/g; #?
		$new_string =~ s/(\w) That (\w)/$1 that $2/g; #?
	}

	$new_string = ucfirst $new_string;

	return ($new_string);

} # End sub make_proper
