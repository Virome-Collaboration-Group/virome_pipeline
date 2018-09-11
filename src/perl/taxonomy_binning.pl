#!/usr/bin/perl

=head1 NAME
   viromeTaxonomyXML.pl

=head1 SYNOPSIS

    USAGE: viromeTaxonomyXML.pl --server server-name --env dbi [--library libraryId]

=head1 OPTIONS

B<--server,-s>
   Server name from where taxonomy info will be collected

B<--library,-l>
    Specific libraryId whoes taxonomy info to collect

B<--env,-e>
    Specific environment where this script is executed.  Based on these values
    db connection and file locations are set.  Possible values are
    igs, dbi, ageek or test

B<--help,-h>
   This help message

=head1  DESCRIPTION
    Create XML document that contains information to draw taxonomy breakdown
    pie chart

=head1  INPUT
    The input is defined with --server,  --library.

=head1  OUTPUT
   XML output per library of taxonomic information

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   viromeTaxonomyXML.pl --server calliope --env dbi --library 31

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
my $cmd = "";

$this->{output_dir} = "$options{outdir}/xDocs";

## make sure everything passed was peachy
&check_parameters(\%options);

#init db connection
my $dbh = DBI->connect("dbi:SQLite:dbname=$options{input}", "", "", { RaiseError => 1}) or die $DBI::errstr;

my $libId = 1;

my $tax_sel = $dbh->prepare(q{SELECT	b.domain,
					b.kingdom,
					b.phylum,
					b.class,
					b.`order`,
					b.family,
					b.genus,
					b.species,
					b.organism,
					b.sequenceId
				FROM	blastp b
				    LEFT JOIN
					sequence s on s.id=b.sequenceId
				WHERE	b.sys_topHit=1
					and b.e_value <= 0.001
					and b.database_name = 'UNIREF100P'
					and s.libraryId = ?
				ORDER BY b.domain,b.kingdom,b.phylum,b.class,b.`order`,b.family,b.genus,b.species});

$tax_sel->execute($libId);
my $rslt = $tax_sel->fetchall_arrayref({});
my $taxStruct;

my $tStruct = {};
foreach my $rec (@$rslt){
	if (!length($rec->{'domain'})){
		$rec->{'domain'} = "UNKNOWN DOMAIN";
	} elsif (!length($rec->{'kingdom'})){
		$rec->{'kingdom'} = "UNKNOWN KINGDOM";
	} elsif (!length($rec->{'phylum'})){
		$rec->{'phylum'} = "UNKNOWN PHYLUM";
	} elsif (!length($rec->{'class'})){
		$rec->{'class'} = "UNKNOWN CLASS";
	} elsif (!length($rec->{'order'})){
		$rec->{'order'} = "UNKNOWN ORDER";
	} elsif (!length($rec->{'family'})){
		$rec->{'family'} = "UNKNOWN FAMILY";
	} elsif (!length($rec->{'genus'})){
		$rec->{'genus'} =  "UNKNOW GENUS";
	} elsif (!length($rec->{'species'})){
		$rec->{'species'} = "UNKNOWN SPECIES";
	} elsif (!length($rec->{'organism'})){
		$rec->{'organism'} = "UNKNOWN ORGANISM";
	}

	$tStruct = createStruct($tStruct, $rec, $rec->{'sequenceId'});
}

my $xml_file = "TAXONOMY_XMLDOC_".$libId.".xml";
my $id_file = "TAXONOMY_IDDOC_".$libId.".xml";

my $xml_out = new IO::File(">". "$this->{output_dir}/$xml_file") or $logger->logdie("Could not open file $this->{output_dir}/$xml_file to write");
my $id_out = new IO::File(">". "$this->{output_dir}/$id_file") or $logger->logdie("Could not open file $this->{output_dir}/$id_file to write");

my $xml_writer = new XML::Writer(OUTPUT=>$xml_out);
my $id_writer = new XML::Writer(OUTPUT=>$id_out);

$xml_writer->xmlDecl("UTF-8");
$id_writer->xmlDecl("UTF-8");

$xml_writer->startTag("root");
$id_writer->startTag("root");

($xml_writer, $id_writer) = printStruct($tStruct, $xml_writer, $id_writer, $id_file);

$xml_writer->endTag("root");
$id_writer->endTag("root");

$xml_writer->end();
$id_writer->end();

$xml_out->close();
$id_out->close();

#### move updated sqlite database back from local to NFS mount
#$cmd = "mv $this->{scratch}/$config->{ToolInfo}->{database_name}->{value} $this->{output_dir}/$config->{ToolInfo}->{database_name}->{value}";
#execute_cmd($cmd);

$logger->info("Taxonomy binning for $options{sample} completed");
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

###############################################################################
sub createStruct{
    my ($taxStruct, $lineage, $seqId) = @_;

    #check for domain. if it exist increment counter by 1.  if not
    #add domain and since the first time, kingdom..3 are also not defined
    #for domain so add that to the struct.
    if (defined $taxStruct->{$lineage->{'domain'}}){
	$taxStruct->{$lineage->{'domain'}}->{'count'}++;
	$taxStruct->{$lineage->{'domain'}}->{'idList'} .= ", $seqId";

	#check for kingdom.  Same logic applies to kingdom
	if (defined $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}){
	    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'count'}++;
	    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'idList'} .= ", $seqId";

	    #check for phylum.  Same logic applies to phylum
	    if (defined $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}){
		$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'count'}++;
		$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'idList'} .= ", $seqId";

		#check for class.
		if (defined $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}){
		    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'count'}++;
		    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'idList'} .= ", $seqId";

		    #check for order
		    if (defined $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}){
			$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'count'}++;
			$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'idList'} .= ", $seqId";

			#check for family
			if (defined $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}){
			    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'count'}++;
			    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'idList'} .= ", $seqId";

			    #check for genus
			    if (defined $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}){
				$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'count'}++;
				$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'idList'} .= ", $seqId";

				#check for species
				if (defined $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'species'}->{$lineage->{'species'}}){
				    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'species'}->{$lineage->{'species'}}->{'count'}++;
				    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'species'}->{$lineage->{'species'}}->{'idList'} .= ", $seqId";

				    #check for organism
				    if (defined $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'species'}->{$lineage->{'species'}}->{'organism'}->{$lineage->{'organism'}}){
					$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'species'}->{$lineage->{'species'}}->{'organism'}->{$lineage->{'organism'}}->{'count'}++;
					$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'species'}->{$lineage->{'species'}}->{'organism'}->{$lineage->{'organism'}}->{'idList'} .= ", $seqId";
				    } else {
    					$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'species'}->{$lineage->{'species'}}->{'organism'}->{$lineage->{'organism'}} = {
																																									'count' => 1,
																																									'idList' => $seqId
																																									}
				    }
				} else {
				    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}}->{'species'}->{$lineage->{'species'}} = {
																																			'count' => 1,
																																			'idList' => $seqId,
																																			'organism' => { $lineage->{'organism'} => {
																																						   'count' => 1,
																																						   'idList' => $seqId
																																						   }
																																			   }
																																		    }
				}
			    } else {
				$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}}->{'genus'}->{$lineage->{'genus'}} = {
																											'count' => 1,
																											'idList' => $seqId,
																									    		'species' => { $lineage->{'species'} => {
																							    								    'count' => 1,
																					    										    'idList' => $seqId,
																			    												    'organism' => { $lineage->{'organism'} => {
																																				   'count' => 1,
															    																					   'idList' => $seqId
													    																							   }
																																	    }
																															    }
																														    }
																											};
			    }
			} else {
			    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}}->{'family'}->{$lineage->{'family'}} = {
																											'count' => 1,
																											'idList' => $seqId,
																											'genus' => { $lineage->{'genus'} => {
																													'count' => 1,
																													'idList' => $seqId,
																													'species' => { $lineage->{'species'} => {
																																	    'count' => 1,
																																	    'idList' => $seqId,
																																	    'organism' => { $lineage->{'organism'} => {
																																						   'count' => 1,
																																						   'idList' => $seqId
																																						   }
																																			   }
																																	    }
																														    }
																													}
																										    }
																								       };
			}
		    } else {
			$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}}->{'order'}->{$lineage->{'order'}} = { 'count' => 1,
																					'idList' => $seqId,
																					 'family' => { $lineage->{'family'} => {
																									'count' => 1,
																									'idList' => $seqId,
																									'genus' => { $lineage->{'genus'} => {
																													'count' => 1,
																													'idList' => $seqId,
																													'species' => { $lineage->{'species'} => {
																																	    'count' => 1,
																																	    'idList' => $seqId,
																																	    'organism' => { $lineage->{'organism'} => {
																																						   'count' => 1,
																																						   'idList' => $seqId
																																						   }
																																			   }
																																	    }
																														    }
																													}
																										    }
																								       }
																						   }
																				        };
		    }
		} else {
		    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}}->{'class'}->{$lineage->{'class'}} = { 'count' => 1,
																	'idList' => $seqId,
																	 'order' => { $lineage->{'order'} => { 'count' => 1,
																					'idList' => $seqId,
																					 'family' => { $lineage->{'family'} => {
																									'count' => 1,
																									'idList' => $seqId,
																									'genus' => { $lineage->{'genus'} => {
																													'count' => 1,
																													'idList' => $seqId,
																													'species' => { $lineage->{'species'} => {
																																	    'count' => 1,
																																	    'idList' => $seqId,
																																	    'organism' => { $lineage->{'organism'} => {
																																						   'count' => 1,
																																						   'idList' => $seqId
																																						   }
																																			   }
																																	    }
																														    }
																													}
																										    }
																								       }
																						   }
																				        }
																		   }
																	}
		}
	    } else {
		$taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}}->{'phylum'}->{$lineage->{'phylum'}} = { 'count' => 1,
													'idList' => $seqId,
													 'class' => { $lineage->{'class'} => { 'count' => 1,
																	'idList' => $seqId,
																	 'order' => { $lineage->{'order'} => { 'count' => 1,
																					'idList' => $seqId,
																					 'family' => { $lineage->{'family'} => {
																									'count' => 1,
																									'idList' => $seqId,
																									'genus' => { $lineage->{'genus'} => {
																													'count' => 1,
																													'idList' => $seqId,
																													'species' => { $lineage->{'species'} => {
																																	    'count' => 1,
																																	    'idList' => $seqId,
																																	    'organism' => { $lineage->{'organism'} => {
																																						   'count' => 1,
																																						   'idList' => $seqId
																																						   }
																																			   }
																																	    }
																														    }
																													}
																										    }
																								       }
																						   }
																				        }
																		   }
																	}
														   }
													};
	    }
	} else {
	    $taxStruct->{$lineage->{'domain'}}->{'kingdom'}->{$lineage->{'kingdom'}} = { 'count' => 1,
									'idList' => $seqId,
									 'phylum' => { $lineage->{'phylum'} => { 'count' => 1,
													'idList' => $seqId,
													 'class' => { $lineage->{'class'} => { 'count' => 1,
																	'idList' => $seqId,
																	 'order' => { $lineage->{'order'} => { 'count' => 1,
																					'idList' => $seqId,
																					 'family' => { $lineage->{'family'} => {
																									'count' => 1,
																									'idList' => $seqId,
																									'genus' => { $lineage->{'genus'} => {
																													'count' => 1,
																													'idList' => $seqId,
																													'species' => { $lineage->{'species'} => {
																																	    'count' => 1,
																																	    'idList' => $seqId,
																																	    'organism' => { $lineage->{'organism'} => {
																																						   'count' => 1,
																																						   'idList' => $seqId
																																						   }
																																			   }
																																	    }
																														    }
																													}
																										    }
																								       }
																						   }
																				        }
																		   }
																	}
														   }
													}
										    }
								       };
	}
    } else {
	$taxStruct->{$lineage->{'domain'}} = { 'count' => 1,
					'idList' => $seqId,
					 'kingdom' => { $lineage->{'kingdom'} => { 'count' => 1,
									'idList' => $seqId,
									'phylum' => { $lineage->{'phylum'} => { 'count' => 1,
													'idList' => $seqId,
													'class' => { $lineage->{'class'} => { 'count' => 1,
																	'idList' => $seqId,
																	'order' => { $lineage->{'order'} => { 'count' => 1,
																					'idList' => $seqId,
																					'family' => { $lineage->{'family'} => {
																									'count' => 1,
																									'idList' => $seqId,
																									'genus' => { $lineage->{'genus'} => {
																													'count' => 1,
																													'idList' => $seqId,
																													'species' => { $lineage->{'species'} => {
																																	    'count' => 1,
																																	    'idList' => $seqId,
																																	    'organism' => { $lineage->{'organism'} => {
																																						   'count' => 1,
																																						   'idList' => $seqId
																																						   }
																																			   }
																																	    }
																														    }
																													}
																										    }
																								       }
																						   }
																				        }
																		   }
																	}
														   }
													}
										    }
								       }
						    }
					};
    }

    return $taxStruct;
}
###############################################################################
sub printStruct{
    my ($tStruct, $xw, $iw, $fname) = @_;

    my $tag = 1;

    for my $domain ( sort keys %$tStruct ) {
	#print "DOMAIN: $domain --> $tStruct->{$domain}->{'count'}\n";
	$iw->emptyTag("TAG_".$tag, 'IDLIST'=>$tStruct->{$domain}->{'idList'});
	$xw->startTag('DOMAIN', 'LABEL'=>$domain,
				    'NAME'=>$domain,
				    'VALUE'=>$tStruct->{$domain}->{'count'},
				    'IDFNAME'=>$fname,
				    'TAG'=>'TAG_'.$tag++);

	my $f2 = $tStruct->{$domain}->{'kingdom'};

	for my $kingdom ( sort keys %$f2 ) {
	    #print "\tFxn2: $kingdom--> $f2->{ $kingdom}->{'count'}\n";
	    $iw->emptyTag("TAG_".$tag, 'IDLIST'=>$f2->{$kingdom}->{'idList'});
	    $xw->startTag('KINGDOM', 'LABEL'=>$kingdom,
					'NAME'=>$kingdom,
					'VALUE'=>$f2->{$kingdom}->{'count'},
					'IDFNAME'=>$fname,
					'TAG'=>'TAG_'.$tag++);

	    my $f3 = $f2->{$kingdom}->{'phylum'};

	    for my $phylum ( sort keys %$f3 ) {
		#print "\t\tFxn3: $phylum --> $f3->{$phylum}->{'count'}\n";
		$iw->emptyTag("TAG_".$tag, 'IDLIST'=>$f3->{$phylum}->{'idList'});
		$xw->startTag('PHYLUM', 'LABEL'=>$phylum,
					    'NAME'=>$phylum,
					    'VALUE'=>$f3->{$phylum}->{'count'},
					    'IDFNAME'=>$fname,
					    'TAG'=>'TAG_'.$tag++);

		my $f4 = $f3->{$phylum}->{'class'};

		for my $class ( sort keys %$f4 ) {
		    #print "\t\t\tFxn4: $class --> $f4->{$class}->{'count'}\n";
		    $iw->emptyTag("TAG_".$tag, 'IDLIST'=>$f4->{$class}->{'idList'});
		    $xw->startTag('CLASS', 'LABEL'=>$class,
						'NAME'=>$class,
						'VALUE'=>$f4->{$class}->{'count'},
						'IDFNAME'=>$fname,
						'TAG'=>'TAG_'.$tag++);

		    my $f5 = $f4->{$class}->{'order'};

		    for my $order ( sort keys %$f5 ) {
			#print "\t\t\t\tFxn5: $order --> $f5->{$order}->{'count'}\n";
			$iw->emptyTag("TAG_".$tag, 'IDLIST'=>$f5->{$order}->{'idList'});
			$xw->startTag('ORDER', 'LABEL'=>$order,
						    'NAME'=>$order,
						    'VALUE'=>$f5->{$order}->{'count'},
						    'IDFNAME'=>$fname,
						    'TAG'=>'TAG_'.$tag++);

			my $f6 = $f5->{$order}->{'family'};

			for my $family ( sort keys %$f6 ) {
			    #print "\t\t\t\t\tFxn6: $family --> $f6->{$family}->{'count'}\n";
			    $iw->emptyTag("TAG_".$tag, 'IDLIST'=>$f6->{$family}->{'idList'});
			    $xw->startTag('FAMILY', 'LABEL'=>$family,
							'NAME'=>$family,
							'VALUE'=>$f6->{$family}->{'count'},
							'IDFNAME'=>$fname,
							'TAG'=>'TAG_'.$tag++);

			    my $f7 = $f6->{$family}->{'genus'};

			    for my $genus ( sort keys %$f7 ) {
				$iw->emptyTag("TAG_".$tag, 'IDLIST'=>$f7->{$genus}->{'idList'});
				$xw->startTag('GENUS', 'LABEL'=>$genus,
							    'NAME'=>$genus,
							    'VALUE'=>$f7->{$genus}->{'count'},
							    'IDFNAME'=>$fname,
							    'TAG'=>'TAG_'.$tag++);

				my $f8 = $f7->{$genus}->{'species'};

				for my $species ( sort keys %$f8 ){
				    $iw->emptyTag("TAG_".$tag, 'IDLIST'=>$f8->{$species}->{'idList'});
				    $xw->startTag('SPECIES', 'LABEL'=>$species,
								'NAME'=>$species,
								'VALUE'=>$f8->{$species}->{'count'},
								'IDFNAME'=>$fname,
								'TAG'=>'TAG_'.$tag++);

				    my $f9 = $f8->{$species}->{'organism'};

				    for my $organism ( sort keys %$f9 ){
					$iw->emptyTag("TAG_".$tag, 'IDLIST'=>$f9->{$organism}->{'idList'});
					$xw->emptyTag('ORGANISM', 'LABEL'=>$organism,
								    'NAME'=>$organism,
								    'VALUE'=>$f9->{$organism}->{'count'},
								    'IDFNAME'=>$fname,
								    'TAG'=>'TAG_'.$tag++);
				    }
				    $xw->endTag('SPECIES');
				}
				$xw->endTag('GENUS');
			    }
			    $xw->endTag('FAMILY');
			}
			$xw->endTag('ORDER');
		    }
		    $xw->endTag('CLASS');
		}
		$xw->endTag('PHYLUM');
	    }
	    $xw->endTag('KINGDOM');
	}
	$xw->endTag('DOMAIN');
    }

    return ($xw, $iw);
}

###############################################################################
