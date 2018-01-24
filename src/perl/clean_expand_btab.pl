#! /usr/bin/perl

=head1 NAME

   expand_uniref100P_btab.pl

=head1 SYNOPSIS

    USAGE: expand_uniref100P_btab.pl --input filename.btab --ouput outputfile.modified.btab

=head1 OPTIONS

B<--input,-i>
   Input file name

B<--help,-h>
   This help message


=head1  DESCRIPTION

    Expand UNIREF100P btab blast output with KEGG, COG, SEED and ACLAME
    results.


=head1  INPUT

    The input is defined with --input.  Input must be a btab blast output
    --output is the full path to output file

    Expected input from RUBBLE in order is

    1 qseqid,qlen,sseqid,salltitles,qstart,qend,sstart,send,pident,ppos,score,bitscore,slen,evalue

    1   query_name (qseqid)
    2   query_length (qlen)
    3   hit_name (sseqid)
    4   hit_description (salltitles)
    5   qry_start (qstart)
    6   qry_end (qend)
    7   hit_start (sstart)
    8   hit_end (send)
    9   percent_identity (pident)
    10  percent_similarity (ppos)
    11  raw_score (score)
    12  bit_score (bitscore)
    13  hit_length (slen)
    14  e_value

=head1  OUTPUT

    Clean btab file, if METAGENOMES hit update hit description

    1   query_name
    2   query_length
    3   algorithm
    4   database_name
    5   db_ranking_code
    6   hit_name
    7   hit_description
    8   qry_start
    9   qry_end
    10  hit_start
    11  hit_end
    12  percent_identity
    13  percent_similarity
    14  raw_score
    15  bit_score
    16  blast frame    # dummy place holder to delete later
    17  query strand   # dummy place holder to delete later
    18  hit_length
    19  e_value

    if UNIREF100P blast expanded btab blast output with
    KEGG, COG, SEED, PHGSEED and ACLAME results and append taxonomy data.

    20  domain
    21  kingdom
    22  phylum
    23  class
    24  order
    25  family
    26  genus
    27  species
    28  organism
	29  functional hit

=head1  CONTACT

  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE

  expand_uniref100P_btab.pl -i input_file_name -ld lookup/file/dir -o output/dir -e igs

=cut

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use UTILS_V;
use DBI;
use Data::Dumper;
use File::Basename;

BEGIN {
  use Ergatis::Logger;
}
##############################################################################
my %options = ();
my $results =
  GetOptions( \%options,
            'input|i=s',
            'subjectDB|s=s',
            'output_dir|o=s',
            'database|b=s',
            'log|l=s',
            'debug|d=s',
            'help|h' )
        || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                 'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

##############################################################################
#### DEFINE GLOBAL VAIRABLES.
##############################################################################

## make sure everything passed was peachy
&check_parameters( \%options );

##############################################################################
my $utils = new UTILS_V;

my $filename = fileparse($options{input}, qr/\.[^.]*/);

open( BTAB, "<", $options{input} ) or $logger->logdie("Can not open file $options{input}");
open( OUT, ">", $options{output_dir} ."/". $filename .".modified.btab" ) or $logger->logdie("Can not open file to write ".$options{output_dir}."/".$filename.".modified.btab");

my $prev         = "";
my $curr         = "";
my @seqarray     = ();
my @unirefarray  = ();
my @keggarray    = ();
my @cogarray     = ();
my @aclamearray  = ();
my @seedarray    = ();
my @phgseedarray = ();

#### database connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{database}", "", "", { RaiseError => 1 }) or $logger->logdie($DBI::errstr);

my $sth = $dbh->prepare('SELECT id,name FROM sequence WHERE typeId=?');

#### if blast against uniref or mgol then sequenceId lookup to orf aa
#### hence typeId=3 else if blast against univec or rna
#### then sequenceId lookup to read ids hence typeId=1

if ($options{subjectDB} =~ /^(uniref100p|metagenomes)$/i) {
    $sth->execute(3);
} else {
    $sth->execute(1);
}

my %sequence_hash;
while ( my $result = $sth->fetchrow_hashref() ) {
	$sequence_hash{$result->{name}} = $result->{id};
}

while (<BTAB>) {
	my $btabline = $_;
	chomp $btabline;

	my @tmp = split( /\t/, $btabline );

    #### insert database name, qry strand and blast frame into the output line
	$btabline = join( "\t", @tmp[0..1], $options{subjectDB}, @tmp[2..11], "0", "0", @tmp[12..$#tmp] );

	if ( $options{subjectDB} =~ /^uniref100p$/i ) {
		$curr = $tmp[0];

		if ( $curr eq $prev ) {
			push( @seqarray, $btabline );
		} else {
			expand();

			#reset array after expansion for new set of seq.
			@seqarray     = ();
			@unirefarray  = ();
			@keggarray    = ();
			@cogarray     = ();
			@seedarray    = ();
			@aclamearray  = ();
			@phgseedarray = ();

			#insert the first new seq info.
			push( @seqarray, $btabline );
			$prev = $curr;
		}    ##END OF ELSE CONDITION
	}
	elsif ( $options{subjectDB} =~ /^metagenomes$/i ) {
		print OUT modifyDescription($btabline) . "\n";
	}
	else {
		print OUT blastn_expand($btabline) . "\n";
	}
}    ##END OF BTAB file

#### expand last set of sequences.
if ( $options{subjectDB} =~ /^UNIREF100P$/i ) {
	expand();
}

close(OUT);

exit(0);

###############################################################################
####  SUBS
###############################################################################
sub check_parameters {
	my $options = shift;

	#### make sure sample_file and output_dir were passed
    my @required = qw(input output_dir subjectDB database);

	foreach my $key (@required) {
		unless ($options{$key}) {
			pod2usage({-exitval => 2,  -message => "ERROR: Input $key not defined", -verbose => 1, -output => \*STDERR});
		    $logger->logdie("No input defined, plesae read perldoc $0\n\n");
        }
	}
}

##############################################################################
sub expand {

    #### each HSP line is expected to be following format (altered from orignal input in main sub routine)
	####
	#### 0:  query_name         : XTG_ctg7180000040098_772_1_1
	#### 1:  query_length       : 257
	#### 2:  database_name      : uniref100p (added in main subroutine)
	#### 3:  hit_name           : fig|6666666.43263.peg.22
	#### 4:  hit_description    : fig|6666666.43263.peg.22 |~|LINEAGE=~DOMAIN:~Viruses;~KINGDOM:~Viruses;~PHYLUM:~Viruses;~CLASS:~Viruses;~ORDER:~Caudovirales;~FAMILY:~Podoviridae;~GENUS:~Podoviridae;~SPECIES:~Puniceispirillum phage HMO-2011;~ORGANISM:~Puniceispirillum phage HMO-2011|~|UNIREF=~ACCESSION:~fig|6666666.43263.peg.22;~DESC:~Protein of Unknown Function;~FXN:~0|~|SEED=~ACCESSION:~fig|6666666.43263.peg.22;~DESC:~UNKNOWN;~FXN:~0|~|PHGSEED=~ACCESSION:~fig|6666666.43263.peg.22;~DESC:~Protein of Unknown Function;~FXN:~0
	#### 5:  qry_start          : 1
	#### 6:  qry_end            : 257
	#### 7:  hit_start          : 1
	#### 8:  hit_end            : 257
	#### 9:  percent_identity   : 57.20
	#### 10: percent_similarity : 77.82
	#### 11: raw_score          : 861
	#### 12: bit_score          : 336
    #### 13: blast frame        : dummy value of blast frame (added in main subroutine)
    #### 14: qry strand         : dummy value of query strand (added in main subroutine)
	#### 15: hit_length         : 300
	#### 16: e_value            : 4e-112

    #### array is a set of hit for a give query sequences.
    #### we should set only one hit as top_hit for each functional database
    #### so init flag outside for each expand (i.e set of query seqeunces)
    my $has_aclame_fxn = 0;
    my $has_cog_fxn = 0;
    my $has_kegg_fxn = 0;
    my $has_seed_fxn = 0;
    my $has_phgseed_fxn = 0;
    my $has_uniref_fxn = 0;

	foreach my $seqline (@seqarray) {
		chomp $seqline;

        my (%data, %taxonomy, %uniref, %seed, %kegg, %cog, %aclame, %phgseed);
        my $lineage = "";

        ####split blast output.
        my @hsp = split(/\t/, $seqline);
        my $sequenceId = $sequence_hash{$hsp[0]};

        #### re-arrange all element as expected by sqlite blastp table
        #### resulting array look like
        ####
        #### 0: sequenceId
        #### 1: qname
        #### 2: qlen
        #### 3: algo
        #### 4: dname
        #### 5: db_ranking
        #### 6: hname
        #### 7: hdesc
        #### 8: qstart
        #### 9: qend
        #### 10: hstart
        #### 11: hend
        #### 12: pident
        #### 13: psim
        #### 14: rscr
        #### 15: bscr
        #### 16: blast frame    # dummy place holder to delete later
        #### 17: query strand   # dummy place holder to delete later
        #### 18: slen
        #### 19: eval
        #### 20: dom
        #### 21: kin
        #### 22: phl
        #### 23: cls
        #### 24: ord
        #### 25: fam
        #### 26: gen
        #### 27: spe
        #### 28: org
        #### 29: fxn_topHit

        my @new = $sequenceId;    # sequenceid add to array
        push (@new, @hsp[0..1]);  # query_name and query_length added
        push (@new, "BLASTP");    # algorithm added
        push (@new, $hsp[2]);     # database name

        #### db ranking system
        if ($hsp[2] =~ /uniref100p/i){
            push (@new, 100);
        } elsif ($hsp[2] =~ /aclame/i){
            push (@new, 90);
        } elsif ($hsp[2] =~ /phgseed/i){
            push (@new, 80);
        } elsif ($hsp[2] =~ /seed/i){
            push (@new, 70);
        } elsif ($hsp[2] =~ /kegg/i){
            push (@new, 60);
        } elsif ($hsp[2] =~ /cog/i){
            push (@new, 50);
        }

        push (@new, @hsp[3..$#hsp]);

        @new = map {$utils->trim($_)} @new;

        #### create a hash with lineage, and each database
        #### resulting hash is of following format
        ####
        #### LINEAGE => DOMAIN:....;KINGDOM:....;PHYLUM:....;......
        #### UNIREF => ACCESSION:.....;DESC:....;FXN:....
        #### SEED => ACCESSION:...;DESC:.....FXN:....
        ####

        my @temp_array = split(/\|~\|/, $hsp[4]);
        @temp_array = @temp_array[1..$#temp_array];
        %data = map_array2hash(\@temp_array, "=~");

        #### extract functional taxonomy
        if (exists $data{LINEAGE}) {
            @temp_array = split(/;~/, $data{LINEAGE});
            %taxonomy = map_array2hash(\@temp_array, ":~");

            $taxonomy{DOMAIN} = "UNKNOWN" unless ($taxonomy{DOMAIN} && $taxonomy{DOMAIN} !~ /null/i);
            $taxonomy{KINGDOM} = "UNKNOWN" unless ($taxonomy{KINGDOM} && $taxonomy{KINGDOM} !~ /null/i);
            $taxonomy{PHYLUM} = "UNKNOWN" unless ($taxonomy{PHYLUM} && $taxonomy{PHYLUM} !~ /null/i);
            $taxonomy{CLASS} = "UNKNOWN" unless ($taxonomy{CLASS} && $taxonomy{CLASS} !~ /null/i);
            $taxonomy{ORDER} = "UNKNOWN" unless ($taxonomy{ORDER} && $taxonomy{ORDER} !~ /null/i);
            $taxonomy{FAMILY} = "UNKNOWN" unless ($taxonomy{FAMILY} && $taxonomy{FAMILY} !~ /null/i);
            $taxonomy{GENUS} = "UNKNOWN" unless ($taxonomy{GENUS} && $taxonomy{GENUS} !~ /null/i);
            $taxonomy{SPECIES} = "UNKNOWN" unless ($taxonomy{SPECIES} && $taxonomy{SPECIES} !~ /null/i);
            $taxonomy{ORGANISM} = "UNKNOWN" unless ($taxonomy{ORGANISM} && $taxonomy{ORGANISM} !~ /null/i);

            $lineage = join("\t", $taxonomy{DOMAIN}, $taxonomy{KINGDOM}, $taxonomy{PHYLUM},
                            $taxonomy{CLASS}, $taxonomy{ORDER}, $taxonomy{FAMILY},
                            $taxonomy{GENUS}, $taxonomy{SPECIES}, $taxonomy{ORGANISM});
        }

        #### expand UNIREF informaiton
        if (exists $data{UNIREF}) {
            #### create a hash of each database value
            #### resulting hash is
            ####
            #### ACCESSION => ....
            #### DESC => ....
            #### FXN => ....

            @temp_array = split(/;~/, $data{UNIREF});
            %uniref = map_array2hash(\@temp_array, ":~");

            $new[4] = "UNIREF100P";
            $new[5] = "100";
            $new[6] = $uniref{ACCESSION};
            $new[7] = ((defined $uniref{DESC}) && (length($uniref{DESC}))) ? $uniref{DESC} : "NA";

            if ($has_uniref_fxn) {
                $uniref{FXN} = 0;
            } else {
                $has_uniref_fxn = $uniref{FXN};
            }

            push (@unirefarray, join("\t", @new, $lineage, $uniref{FXN}));
        }

        #### expand SEED informaiton
        if (exists $data{SEED}) {
            @temp_array = split(/;~/, $data{SEED});
            %seed = map_array2hash(\@temp_array, ":~");

            ###############################################################
            #### a work around for incorrect seed mldbm file when running
            #### uniref_SEPT2013_edit was created in DEC2014.
            ###############################################################
            #my $temp_acc_hash = $lookup_util->get_acc_from_lookup( "seed", $seed{ACCESSION} );
            #my $temp_hash     = $temp_acc_hash->{acc_data}[0];
            #$seed{DESC} = $s_hash->{desc};
            #if ((defined $s_hash->{fxn1}) && length($s_hash->{fxn1})
            #	&& ($s_hash->{fxn1} !~
            #		/unknown|unclassified|unassigned|uncharacterized/i )) {
            #	$seed{FXN} = 1;
            #}

            $new[4] = "SEED";
            $new[5] = "70";
            $new[6] = $seed{ACCESSION};
            $new[7] = ((defined $seed{DESC}) && (length($seed{DESC}))) ? $seed{DESC} : "NA";

            if ($has_seed_fxn) {
                $seed{FXN} = 0;
            } else {
                $has_seed_fxn = $seed{FXN};
            }

            push (@seedarray, join("\t", @new, $lineage, $seed{FXN}));
        }

        #### expand KEGG informaiton
        if (exists $data{KEGG}) {
            @temp_array = split(/;~/, $data{KEGG});
            %kegg = map_array2hash(\@temp_array, ":~");

            $new[4] = "KEGG";
            $new[5] = "60";
            $new[6] = $kegg{ACCESSION};
            $new[7] = ((defined $kegg{DESC}) && (length($kegg{DESC}))) ? $kegg{DESC} : "NA";

            if ($has_kegg_fxn) {
                $kegg{FXN} = 0;
            } else {
                $has_kegg_fxn = $kegg{FXN};
            }

            push (@keggarray, join("\t", @new, $lineage, $kegg{FXN}));
        }

        #### expand COG informaiton
        if (exists $data{COG}) {
            @temp_array = split(/;~/, $data{COG});
            %cog = map_array2hash(\@temp_array, ":~");

            $new[4] = "COG";
            $new[5] = "50";
            $new[6] = $cog{ACCESSION};
            $new[7] = ((defined $cog{DESC}) && (length($cog{DESC}))) ? $cog{DESC} : "NA";

            if ($has_cog_fxn) {
                $cog{FXN} = 0;
            } else {
                $has_cog_fxn = $cog{FXN};
            }

            push (@cogarray, join("\t", @new, $lineage, $cog{FXN}));
        }

        #### expand ACLAME informaiton
        if (exists $data{ACLAME}) {
            @temp_array = split(/;~/, $data{ACLAME});
            %aclame = map_array2hash(\@temp_array, ":~");

            $new[4] = "ACLAME";
            $new[5] = "90";
            $new[6] = $aclame{ACCESSION};
            $new[7] = ((defined $aclame{DESC}) && (length($aclame{DESC}))) ? $aclame{DESC} : "NA";

            if ($has_aclame_fxn) {
                $aclame{FXN} = 0;
            } else {
                $has_aclame_fxn = $aclame{FXN};
            }

            push (@aclamearray, join("\t", @new, $lineage, $aclame{FXN}));
        }

        #### expand PHGSEED informaiton
        if (exists $data{PHGSEED}) {
            @temp_array = split(/;~/, $data{PHGSEED});
            %phgseed = map_array2hash(\@temp_array, ":~");

            $new[4] = "PHGSEED";
            $new[5] = "80";
            $new[6] = $phgseed{ACCESSION};
            $new[7] = ((defined $phgseed{DESC}) && (length($phgseed{DESC}))) ? $phgseed{DESC} : "NA";

            if ($has_phgseed_fxn) {
                $phgseed{FXN} = 0;
            } else {
                $has_phgseed_fxn = $phgseed{FXN};
            }

            push (@phgseedarray, join("\t", @new, $lineage, $phgseed{FXN}));
        }
    }

    ####PRINT UNIREF FIRST
    foreach my $unirefline (@unirefarray) {
        print OUT $unirefline."\n";
    }

    foreach my $aclameline (@aclamearray) {
        print OUT $aclameline."\n";
    }

    foreach my $phgseedline (@phgseedarray) {
        print OUT $phgseedline."\n";
    }

    foreach my $seedline (@seedarray) {
        print OUT $seedline."\n";
    }

    foreach my $keggline (@keggarray) {
        print OUT $keggline."\n";
    }

    foreach my $cogline (@cogarray) {
        print OUT $cogline."\n";
    }
}

##############################################################################
sub modifyDescription {
    my $seqline = shift;

	my @hsp = split(/\t/, $seqline);
    my $sequenceId = $sequence_hash{$hsp[0]};
	my $str = "";

    #### re-arrange all element as expected by sqlite blastp table
    #### resulting array look like
    ####
    #### 0: sequenceId
    #### 1: qname
    #### 2: qlen
    #### 3: algo
    #### 4: dname
    #### 5: db_ranking
    #### 6: hname
    #### 7: hdesc
    #### 8: qstart
    #### 9: qend
    #### 10: hstart
    #### 11: hend
    #### 12: pident
    #### 13: psim
    #### 14: rscr
    #### 15: bscr
    #### 16: blast frame    # dummy place holder to delete later
    #### 17: query strand   # dummy place holder to delete later
    #### 18: slen
    #### 19: eval

    my @new = $sequenceId;    # sequenceid add to array
    push (@new, @hsp[0..1]);  # query_name and query_length added
    push (@new, "BLASTP");    # algorithm added
    push (@new, $hsp[2]);     # database name
    push (@new, 10);          # db ranking system

    push (@new, @hsp[3..$#hsp]);

    @new = map {$utils->trim($_)} @new;

	return join("\t", @new);
}

##############################################################################
sub blastn_expand {
    my $seqline = shift;

	my @hsp = split(/\t/, $seqline);
    my $sequenceId = $sequence_hash{$hsp[0]};
	my $str = "";

    #### re-arrange all element as expected by sqlite blastp table
    #### resulting array look like
    ####
    #### 0: sequenceId
    #### 1: qname
    #### 2: qlen
    #### 3: algo
    #### 4: dname
    #### 5: db_ranking
    #### 6: hname
    #### 7: hdesc
    #### 8: qstart
    #### 9: qend
    #### 10: hstart
    #### 11: hend
    #### 12: pident
    #### 13: psim
    #### 14: rscr
    #### 15: bscr
    #### 16: blast frame    # dummy place holder to delete later
    #### 17: query strand   # dummy place holder to delete later
    #### 18: slen
    #### 19: eval

    my @new = $sequenceId;    # sequenceid add to array
    push (@new, @hsp[0..1]);  # query_name and query_length added
    push (@new, "BLASTN");    # algorithm added
    push (@new, $hsp[2]);     # database name

    if ($hsp[2] =~ /univec/i) {
        push (@new, 4);          # db ranking system
    } else {
        push (@new, 5);          # db ranking system
    }


    push (@new, @hsp[3..$#hsp]);

    @new = map {$utils->trim($_)} @new;

	return join("\t", @new);
}

#############################################################################
sub modify_old_delim {
	my $str = shift;

	$str =~ s/\|\|/\|~\|/g;

	$str =~ s/PHGSEED=ACCESSION:/PHGSEED=~ACCESSION:~/ig;
	$str =~ s/SEED=ACCESSION:/SEED=~ACCESSION:~/ig;
	$str =~ s/KEGG=ACCESSION:/KEGG=~ACCESSION:~/ig;
	$str =~ s/COG=ACCESSION:/COG=~ACCESSION:~/ig;
	$str =~ s/ACLAME=ACCESSION:/ACLAME=~ACCESSION:~/ig;
	$str =~ s/UNIREF=ACCESSION:/UNIREF=~ACCESSION:~/ig;


	$str =~ s/;DESC:/;~DESC:~/ig;
	$str =~ s/;FXN:/;~FXN:~/ig;

	$str =~ s/LINEAGE=DOMAIN:/LINEAGE=~DOMAIN:~/ig;
	$str =~ s/;KINGDOM:/;~KINGDOM:~/ig;
	$str =~ s/;PHYLUM:/;~PHYLUM:~/ig;
	$str =~ s/;CLASS:/;~CLASS:~/ig;
	$str =~ s/;ORDER:/;~ORDER:~/ig;
	$str =~ s/;FAMILY:/;~FAMILY:~/ig;
	$str =~ s/;GENUS:/;~GENUS:~/ig;
	$str =~ s/;SPECIES:/;~SPECIES:~/ig;
	$str =~ s/;ORGANISM:/;~ORGANISM:~/ig;

	return $str;
}

#############################################################################
sub map_array2hash {
	my $array = shift;
	my $delim = shift;

	my @n = @{$array};

	my %hash;

	foreach my $elm (@n) {
		my @a = split(/$delim/, $elm);

		$hash{$a[0]} = $a[1];
	}

	return %hash;
}
