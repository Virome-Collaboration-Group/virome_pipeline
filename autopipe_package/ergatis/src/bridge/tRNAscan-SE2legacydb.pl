#!/usr/local/bin/perl

=head1  NAME 

tRNAscan-SE2legacydb.pl - load tRNAscan-SE data into the legacy db

=head1 SYNOPSIS

USAGE: tRNAscan-SE2legacydb.pl 
        --input_list=/path/to/somefile.raw.list
        --database=aa1
      [ --log=/path/to/some.log 
        --delete_existing=1
      ]

=head1 OPTIONS

B<--input_list,-i> 
    raw list file from a tRNAscan-SE workflow component run.
	
B<--database,-d> 
    Sybase project database ID.

B<--delete_existing,-x> 
    optional.  will first delete all existing tRNA features for each assembly passed.

B<--log,-d> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

This script is used to load tRNAscan-SE evidence into the legacy database schema.
It supports input in raw tRNAscan-SE output format only.

=head1 INPUT

The file names within the raw list must be contain the portion 
'assembly.$assemblyid.$program.raw', like:

    /some/path/aa1.assembly.24963.tRNAscan-SE.raw
    
So that the assembly number (24963) and program can be extracted from the name.

=head1 OUTPUT

This is a database loading script.  There is no other output unless you use 
the --log option.

=head1 CONTACT

    Brett Whitty
	bwhitty@tigr.org

=cut

use strict;
use warnings;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
			  'input_list|i=s',
              'database|d=s',
              'delete_existing|x',
              'log|l=s',
			  'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
check_parameters(\%options);

## get the current user
my ($user) = getpwuid($<);

## a global counter used for numbering tRNA exons
my $NUMBER = 0;

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## connect to the database
my $dbh = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092", 'egc', 'egcpwd');
$dbh->do("use $options{database}");

## load the file list
my @files;
open(my $listfh, "<$options{input_list}") || die "can't open file list: $!";

while (<$listfh>) {
	chomp;
    next if (/^\s*$/);
    _log("adding $_ to the process list");
    push @files, $_;
}
close $listfh;

foreach my $file (@files) {
    _log("processing $file");
    
    my ($asmbl_id, $prog_name);
    if ($file =~ /\.assembly\.(\d+)\.(.+)\.raw/) {
        ($asmbl_id, $prog_name) = ($1, $2);
        _log("extracted asmbl_id $asmbl_id, prog_name $prog_name from $file");
    } else {
        die "unable to parse asmbl_id from $file\n";
    }
	
	## parse the raw file
	_log("parsing $file");
	my $results_array_ref = parse_trnascan_raw($file);
   
	## delete existing tRNA info for the assembly
	if ($options{'delete_existing'} ) {
		do_deletes($dbh, $asmbl_id);
	}
	
	my $ntseq = get_seq($dbh, $asmbl_id);
	
	foreach my $result_line(@{$results_array_ref}) {
		if ($result_line->[4] ne 'Pseudo') { ## skip Pseudo-tRNAs
   			insert_tRNA(	$dbh,
							$result_line->[2],	#start pos
							$result_line->[3],	#end pos
							$result_line->[4],  #aa
							$result_line->[5],  #anti-codon
							$result_line->[6],  #bound1	
							$result_line->[7],	#bound2
							$result_line->[8],	#score
							$asmbl_id,			#assembly id
							$ntseq,				#nt sequence of assembly
						);
		}
	}

}
$dbh->disconnect();

exit();

sub parse_trnascan_raw {
	my $file = shift @_;
	my @results;
	open (IN, $file) || die "can't open raw file '$file' for reading";
	while (<IN>) {
		chomp;
		push(@results, [split("\t")]);
	}
	return \@results;
}

sub _log {
    my $msg = shift;
    
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    
    ## database and input_list are required
    unless ( defined $options{database} && $options{input_list} ) { 
        print STDERR "database and input_list options are required\n\n";
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }
    
    ## make sure input list exists
    unless ( -e $options{input_list} ) {
        print STDERR "\n\ninput_list $options{input_list} does not exist\n\n";
        exit(1);
    }
    
    ## set some defaults
    $options{delete_existing} = 0 unless ( $options{delete_existing} );
}

sub do_deletes {
		my($dbh,$asmbl_id) = @_;
		delete_ALL_features($dbh, $asmbl_id, "tRNA");
		delete_ALL_features($dbh, $asmbl_id, "pre-tRNA");
		delete_ALL_features($dbh, $asmbl_id, "rna-exon");
		$dbh->do("delete feat_link where child_feat like \"$asmbl_id.x%\"");
}

sub delete_ALL_features {  ## Be careful... this removes EVERYTHING!!
                           ## use really for tRNA and ORF
   my($dbh, $asmbl, $type) = @_;
   my($query);
   my($db_query);
  
	_log("Doing '$type' deletes for '$asmbl'...");
   
   $dbh->do("delete evidence from asm_feature f, evidence e"
   		  . " where f.feat_name=e.feat_name"
   		  . " and f.asmbl_id = $asmbl"
   		  . " and f.feat_type = \"$type\"");
   
   $dbh->do("delete ORF_attribute from asm_feature f, ORF_attribute o"
   		  . " where f.feat_name=o.feat_name"
   		  . " and f.asmbl_id = $asmbl"
   		  . " and f.feat_type = \"$type\"");

   $dbh->do("delete phys_ev from asm_feature f, phys_ev o"
   		  . " where f.feat_name=o.feat_name"
   		  . " and f.asmbl_id = $asmbl"
   		  . " and f.feat_type = \"$type\"");
   
   $dbh->do("delete ident from asm_feature f, ident i"
   		  . " where f.feat_name=i.feat_name"
   		  . " and f.asmbl_id = $asmbl"
   		  . " and f.feat_type = \"$type\"");
   
   $dbh->do("delete role_link from role_link r, asm_feature a"
   		  . " where a.asmbl_id = $asmbl"
   		  . " and a.feat_type = \"$type\""
   		  . " and a.feat_name = r.feat_name");
   
   $dbh->do("delete feat_link from feat_link r, asm_feature a"
   		  . " where a.asmbl_id = $asmbl"
   		  . " and a.feat_type = \"$type\""
   		  . " and a.feat_name = r.child_feat");
   
   $dbh->do("delete feat_link from feat_link r, asm_feature a"
   		  . " where a.asmbl_id = $asmbl"
   		  . " and a.feat_type = \"$type\""
   		  . " and a.feat_name = r.parent_feat");
   
   $dbh->do("delete asm_feature from asm_feature"
	   . " where asm_feature.asmbl_id = $asmbl"
       . " and asm_feature.feat_type = \"$type\"");
   
}

sub insert_tRNA {
    my ($dbh, $end5, $end3, $aa, $anti_codon, $bound1, $bound2, $score, $asmbl, $ntseq) = @_;
    my ($att_id, $ret, $prog, @check, $count, $id, $preseq, $seq, @names, $file, $feat_name);
    my ($zerobuf, $i);
    my ($db_query);
    
	my $zeros = "00000";
   	
    $feat_name = "tRNA-" . $aa;

    _log("INSERTING $feat_name...");
    _log("$end5, $end3, $bound1, $bound2");
    ###### create unique name by checking the molecule 
    ###### (would like to do this project wide, but can't right) now for
    ###### the existance of that tRNA already.  Count the 
    ###### number of occurances and increase the count by
    ###### one to generate the new name.
 
    my $query = "select COUNT(feat_name) from asm_feature"
			  . " where feat_name like ?"
	    	  . " and asmbl_id = ?";

	$db_query = $dbh->prepare($query);
	$db_query->execute("$asmbl.$feat_name%", $asmbl);
	my $row = $db_query->fetchrow_arrayref();
	$count = $row->[0];
	$db_query->finish();
 
    $count++;
	
    my $feat_name1 = $asmbl.".pre-".$feat_name . "-" . $count;
    my $feat_name2 = $asmbl.".".$feat_name . "-" . $count;

    _log("## pre-tRNA sequence....");
    $preseq = &cut_seq2($ntseq,$end5,$end3);
    my ($exon1_3,$exon2_5)='';
    my ($e1seq,$e2seq)='';
    if ($bound1 != 0) {
		if ($end5 < $end3) {
	    	$exon1_3 = $bound1 - 1;
		    $exon2_5 = $bound2+1;
		} else {
	    	$exon1_3 = $bound1 + 1;
		    $exon2_5 = $bound2 - 1;
		}
		_log("## tRNA sequence....");
		$e1seq = cut_seq2($ntseq,$end5,$exon1_3);
        $e2seq = cut_seq2($ntseq,$exon2_5,$end3);
		$seq = $e1seq . $e2seq;
    } else {
		_log("## tRNA sequence....");
		$seq = $preseq;
    }
    
	###### load the feature as type pre-tRNA
    $dbh->do("insert asm_feature(feat_name, end5, end3, feat_type, assignby, date, asmbl_id, change_log, save_history, sequence)"
		   . " values (\"$feat_name1\", $end5, $end3, \"pre-tRNA\", \"$user\", getdate(), $asmbl, 0, 1, \"$preseq\")");
    
    ###### load the feature as type tRNA (intron omitted)
    $dbh->do("insert asm_feature(feat_name, end5, end3, feat_type, assignby, date, asmbl_id, change_log, save_history, sequence)"
		   . " values (\"$feat_name2\", $end5, $end3, \"tRNA\", \"$user\", getdate(), $asmbl, 0, 1, \"$seq\")");
	
    ### link the pre-tRNA to the tRNA.
    insert_feat_link($dbh, $feat_name2, $feat_name1);
    
    ##### load the exons
    $NUMBER++;
	
	$zerobuf = substr($zeros, 0, length($zeros) - length($NUMBER));
    my $feat_name3 = $asmbl . ".x". $zerobuf . $NUMBER;
    
    if($bound1 == 0) {
    	$dbh->do("insert asm_feature(feat_name, end5, end3, feat_type, assignby, date, asmbl_id, change_log, save_history, sequence)"
		   	   . " values (\"$feat_name3\", $end5, $end3, \"rna-exon\", \"$user\", getdate(), $asmbl, 0, 1, \"$seq\")");
	
		insert_feat_link($dbh, $feat_name3, $feat_name2);
    } else {
		$NUMBER++;
		
		$zerobuf = substr($zeros, 0, length($zeros) - length($NUMBER));
		my $feat_name4 = $asmbl . ".x". $zerobuf . $NUMBER;
	
    	$dbh->do("insert asm_feature(feat_name, end5, end3, feat_type, assignby, date, asmbl_id, change_log, save_history, sequence)"
		   	   . " values (\"$feat_name3\", $end5, $exon1_3, \"rna-exon\", \"$user\", getdate(), $asmbl, 0, 1, \"$e1seq\")");
		
		insert_feat_link($dbh, $feat_name3,$feat_name2);

    	$dbh->do("insert asm_feature(feat_name, end5, end3, feat_type, assignby, date, asmbl_id, change_log, save_history, sequence)"
		   	   . " values (\"$feat_name4\", $exon2_5, $end3, \"rna-exon\", \"$user\", getdate(), $asmbl, 0, 1, \"$e2seq\")");
	
		insert_feat_link($dbh, $feat_name4, $feat_name2);
    }
    
#	$db_query->finish();
	
    ###### load the ORF_attribute (link by the "model" tRNA and not pre-tRNA)
    $att_id = insert_ORF_attribute($dbh, $feat_name2, "tRNA", 1, "workflow", "$user");
    
    ###### load the scores.
    $anti_codon =~ s/\n//g;
    insert_ORFattribute_tRNA_score ($dbh, $att_id, $anti_codon, $score);

    ## insert ident
    $dbh->do("insert ident (feat_name, com_name) values (\"$feat_name2\", \"$feat_name (anticodon: $anti_codon)\")");
}

sub cut_seq2 {
    my($sequence, $end5, $end3) = @_;
    my($s, $seq_len, $query, $x, $ret, $tmp);
    
	if ($end5 < $end3) {
	  	$s=substr($sequence, $end5 - 1, $end3 - $end5 + 1);
    } else {
  		$tmp=substr($sequence, $end3 - 1, $end5 - $end3 + 1);
	  	$s = join("",reverse(split(/ */,$tmp)));
  		$s=~tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/;
    }
    return($s);
}

sub insert_feat_link {
    my($dbh,$child,$parent) = @_;
    my($query, $db_query);
    
	$dbh->do("insert feat_link (parent_feat, child_feat, assignby, datestamp) "
    	   . "values (\"$parent\", \"$child\", \"$user\", getdate())");
   
}

sub insert_ORF_attribute {
    my($dbh, $feat_name, $type, $curated, $method) = @_;
    my($query, $db_query);

    $dbh->do("insert ORF_attribute (feat_name, att_type, curated, method, date, assignby) "
		   . "values (\"$feat_name\",\"$type\",$curated,\"$method\",getdate(),\"$user\")");
    
	$query = "select id from ORF_attribute "
		   . "where feat_name = ? "
		   . "and att_type = ?";
    
	$db_query = $dbh->prepare($query);
	$db_query->execute($feat_name, $type);
	my $row = $db_query->fetchrow_arrayref();
	my $id = $row->[0];
	$db_query->finish();
	
	return $id;
}

sub insert_ORFattribute_tRNA_score {
    my ($dbh, $id, $anti_codon, $score) = @_;
    my($query, $db_query);
    $dbh->do("update ORF_attribute set "
		   . "score = \'$anti_codon\', "
	       . "score_desc = \'anti-codon\', "
	       . "score2 = \'$score\', "
	       . "score2_desc = \'cove\' "
	       . "where id = $id");
}

sub get_seq {
    my($dbh, $id) = @_;
    my($query, $db_query, $seq_len);
    
    ###### query 1
    $query = "select datalength(sequence) from assembly"
           . " where asmbl_id = $id";
	
	$db_query = $dbh->prepare($query);
	$db_query->execute();
	my $row = $db_query->fetchrow_arrayref();
	$seq_len = $row->[0];
	$db_query->finish();
    
    $seq_len = ($seq_len < 1e06) ? 1e06 : $seq_len;
    
    ###### query 2
    if($seq_len ne ""){
		$dbh->do("set textsize $seq_len");    
    } else {
		_log("ERROR:\tThere is no length associated with the following ASMBL_ID:\t$id");
    }
    ###### query 3
    $query = "select sequence from assembly"
    	   . " where asmbl_id=$id";    
	$db_query = $dbh->prepare($query);
	$db_query->execute();
	
	$row = $db_query->fetchrow_arrayref();
	my $seq = $row->[0];
	$db_query->finish();
    
    $seq =~ s/\s//g;
    $seq = uc $seq;
    return($seq);
}
